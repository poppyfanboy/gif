/**
 *  @typedef {('srgb' | 'linear' | 'lab' | 'oklab')} ColorSpace
 *  @typedef {('none' | 'floyd-steinberg' | 'ordered')} DitheringAlgorithm
 *  @typedef {('fixed' | 'generated')} PaletteType
 *  @typedef {('web-safe' | 'monochrome' | 'black-white')} FixedPalette
 *  @typedef {('median-cut' | 'k-means' | 'k-means++' | 'modified-median-cut' | 'octree')} PaletteAlgorithm
 */

const BYTES_PER_PAGE = 64 * 1024;
const INITIAL_PAGE_COUNT = 16;
const MAX_PAGE_COUNT = 8192;

const ALIGNMENT = 16;

class Arena {
    constructor(wasmMemory, heapBase) {
        this.currentPageCount = INITIAL_PAGE_COUNT;

        this.wasmMemory = wasmMemory;
        this.memory = new DataView(this.wasmMemory.buffer);

        const arenaBase = heapBase + 2 * 4;

        this.beginPointer = heapBase;
        this.memory.setInt32(this.beginPointer, arenaBase, true);

        this.endPointer = heapBase + 4;
        this.memory.setInt32(this.endPointer, this.currentPageCount * BYTES_PER_PAGE, true);
    }

    get begin() {
        return this.memory.getInt32(this.beginPointer, true);
    }

    get end() {
        return this.memory.getInt32(this.endPointer, true);
    }

    get pointer() {
        return this.beginPointer;
    }

    snapshot() {
        return this.begin;
    }

    rewind(snapshot) {
        this.memory.setInt32(this.beginPointer, snapshot, true);
    }

    ensure(size) {
        const memoryLeft = this.end - this.begin;
        if (size <= memoryLeft) {
            return;
        }

        const pagesLeft = MAX_PAGE_COUNT - this.currentPageCount;
        if (memoryLeft + pagesLeft * BYTES_PER_PAGE < size) {
            throw new Error('Out of memory');
        }

        const pagesToCommit = Math.min(pagesLeft, Math.ceil(size / BYTES_PER_PAGE));

        this.wasmMemory.grow(pagesToCommit);
        this.currentPageCount += pagesToCommit;

        this.memory = new DataView(this.wasmMemory.buffer);
        this.memory.setInt32(this.endPointer, this.currentPageCount * BYTES_PER_PAGE, true);
    }

    alloc(size) {
        if (size == 0) {
            return 0;
        }

        const padding = (~this.begin + 1) & (ALIGNMENT - 1);
        this.ensure(padding + size);

        const newPointer = this.begin + padding;
        this.memory.setInt32(this.beginPointer, this.begin + padding + size, true);
        return newPointer;
    }
}

const wasmMemory = new WebAssembly.Memory({
    initial: INITIAL_PAGE_COUNT,
    maximum: MAX_PAGE_COUNT,
});

const gifEncoderPromise = (async () => {
    const wasmModule = await WebAssembly.instantiateStreaming(fetch('gif_encoder.wasm'), {
        env: {
            memory: wasmMemory,
            reportAlloc: (_arenaPointer, size) => {
                const padding = (~arena.begin + 1) & (ALIGNMENT - 1);
                arena.ensure(padding + size);
            },

            pow: (base, exponent) => Math.pow(base, exponent),
            cbrt: (value) => Math.cbrt(value),
            round: (value) => Math.round(value),
        },
    });

    const arena = new Arena(wasmMemory, wasmModule.instance.exports.get_heap_base());

    const {
        srgb_palette_black_and_white,
        srgb_palette_monochrome,
        srgb_palette_web_safe,

        srgb_to_float, srgb_to_linear, srgb_to_lab, srgb_to_oklab,
        float_to_srgb, linear_to_srgb, lab_to_srgb, oklab_to_srgb,

        colors_unique_inplace,

        palette_by_median_cut,
        palette_by_modified_median_cut,
        palette_by_k_means,
        palette_by_octree,

        image_quantize,
        image_floyd_steinberg_dither,
        image_ordered_dither,

        gif_out_buffer_create,
        gif_out_buffer_capacity_left,
        gif_out_buffer_grow,

        gif_encoder_create,
        gif_encoder_start,
        gif_encoder_finish,

        gif_encoder_start_frame,
        gif_encoder_feed_frame,
        gif_encoder_finish_frame,
    } = wasmModule.instance.exports;

    /**
    * @param {number} srgbColors Pointer to the array of bytes.
    * @param {number} colorCount
    * @param {number} components
    * @param {ColorSpace} colorSpace
    *
    * @returns {number} Pointer to the arena-allocated color array.
    */
    function fromSrgb(srgbColors, colorCount, components, colorSpace) {
        switch (colorSpace) {
        case 'srgb':    return srgb_to_float(srgbColors, colorCount, components, arena.pointer);
        case 'linear':  return srgb_to_linear(srgbColors, colorCount, components, arena.pointer);
        case 'lab':     return srgb_to_lab(srgbColors, colorCount, components, arena.pointer);
        case 'oklab':   return srgb_to_oklab(srgbColors, colorCount, components, arena.pointer);
        }
    }

    /**
    * @param {number} colors Pointer to the array of floats.
    * @param {number} colorCount
    * @param {number} components
    * @param {ColorSpace} colorSpace
    *
    * @returns {number} Pointer to the arena-allocated color array.
    */
    function intoSrgb(colors, colorCount, components, colorSpace) {
        switch (colorSpace) {
        case 'srgb':    return float_to_srgb(colors, colorCount, components, arena.pointer);
        case 'linear':  return linear_to_srgb(colors, colorCount, components, arena.pointer);
        case 'lab':     return lab_to_srgb(colors, colorCount, components, arena.pointer);
        case 'oklab':   return oklab_to_srgb(colors, colorCount, components, arena.pointer);
        }
    }

    /**
    * @param {GifConfig} config
    */
    function encode(image, config) {
        const arenaSnapshot = arena.snapshot();

        const srgbImage = arena.alloc(image.width * image.height * image.components);
        {
            // Don't keep typed arrays in your hands for too long: as the memory grows they rot
            // away.
            const srgbImageArray = new Uint8Array(
                arena.memory.buffer,
                srgbImage, image.width * image.height * image.components
            );
            srgbImageArray.set(image.data, 0);
        }

        const floatImage = fromSrgb(
            srgbImage, image.width * image.height, image.components,
            config.colorSpace
        );

        let colorCount;
        const colorCountPointer = arena.alloc(4);

        let srgbColors, floatColors;

        if (config.paletteType == 'fixed') {
            const colorCountPointer = arena.alloc(4);

            if (config.fixedPalette == 'web-safe') {
                srgbColors = srgb_palette_web_safe(
                    colorCountPointer, image.components, arena.pointer
                );
            } else if (config.fixedPalette == 'monochrome') {
                srgbColors = srgb_palette_monochrome(
                    colorCountPointer, image.components, arena.pointer
                );
            } else if (config.fixedPalette == 'black-white') {
                srgbColors = srgb_palette_black_and_white(
                    colorCountPointer, image.components, arena.pointer
                );
            }

            colorCount = arena.memory.getInt32(colorCountPointer, true);
            floatColors = fromSrgb(srgbColors, colorCount, image.components, config.colorSpace);
        } else {
            const uniquePixelCountPointer = arena.alloc(4);
            const uniqueSrgbPixels = colors_unique_inplace(
                srgbImage, image.width * image.height, image.components,
                uniquePixelCountPointer,
                arena.pointer
            );
            const uniquePixelCount = arena.memory.getInt32(uniquePixelCountPointer, true);

            const uniquePixels = fromSrgb(
                uniqueSrgbPixels, uniquePixelCount, image.components,
                config.colorSpace
            );


            switch (config.paletteAlgorithm) {
            case 'median-cut': {
                floatColors = palette_by_median_cut(
                    uniquePixels, uniquePixelCount, image.components,
                    config.paletteColorCount,
                    colorCountPointer,
                    arena.pointer
                );
                colorCount = arena.memory.getInt32(colorCountPointer, true);

                srgbColors = intoSrgb(floatColors, colorCount, image.components, config.colorSpace);
            } break;

            case 'k-means': {
                floatColors = palette_by_k_means(
                    uniquePixels, uniquePixelCount, image.components,
                    config.paletteColorCount,
                    colorCountPointer,
                    false,
                    arena.pointer
                );
                colorCount = arena.memory.getInt32(colorCountPointer, true);

                srgbColors = intoSrgb(floatColors, colorCount, image.components, config.colorSpace);
            } break;

            case 'k-means++': {
                floatColors = palette_by_k_means(
                    floatImage, image.width * image.height, image.components,
                    config.paletteColorCount,
                    colorCountPointer,
                    true,
                    arena.pointer
                );
                colorCount = arena.memory.getInt32(colorCountPointer, true);

                srgbColors = intoSrgb(floatColors, colorCount, image.components, config.colorSpace);
            } break;

            case 'modified-median-cut': {
                srgbColors = palette_by_modified_median_cut(
                    uniqueSrgbPixels, uniquePixelCount, image.components,
                    config.paletteColorCount,
                    colorCountPointer,
                    arena.pointer
                );
                colorCount = arena.memory.getInt32(colorCountPointer, true);

                floatColors = fromSrgb(srgbColors, colorCount, image.components, config.colorSpace);
            } break;

            case 'octree': {
                srgbColors = palette_by_octree(
                    uniqueSrgbPixels, uniquePixelCount, image.components,
                    config.paletteColorCount,
                    colorCountPointer,
                    arena.pointer
                );
                colorCount = arena.memory.getInt32(colorCountPointer, true);

                floatColors = fromSrgb(srgbColors, colorCount, image.components, config.colorSpace);
            } break;
            }
        }

        let indexedImage;
        if (config.ditheringAlgorithm == 'floyd-steinberg') {
            indexedImage = image_floyd_steinberg_dither(
                floatImage, image.width, image.height, image.components,
                floatColors, colorCount,
                arena.pointer
            );
        } else if (config.ditheringAlgorithm == 'ordered') {
            indexedImage = image_ordered_dither(
                floatImage, image.width, image.height, image.components,
                floatColors, colorCount,
                arena.pointer
            );
        } else {
            indexedImage = image_quantize(
                floatImage, image.width * image.height, image.components,
                floatColors, colorCount,
                arena.pointer
            );
        }

        let outputBufferCapacity = 64 * 1024;
        const outputBuffer = gif_out_buffer_create(outputBufferCapacity, arena.pointer);
        const gifEncoder = gif_encoder_create(arena.pointer);

        gif_encoder_start(
            gifEncoder,
            image.width, image.height,
            srgbColors, colorCount, image.components,
            outputBuffer
        );

        let indexedImageIter = indexedImage;
        const indexedImageEnd = indexedImage + image.width * image.height;

        gif_encoder_start_frame(gifEncoder, 0, 0, image.components, outputBuffer);

        while (indexedImageIter != indexedImageEnd) {
            indexedImageIter += gif_encoder_feed_frame(
                gifEncoder,
                indexedImageIter,
                indexedImageEnd - indexedImageIter,
                outputBuffer
            );
            if (gif_out_buffer_capacity_left(outputBuffer) < 1024) {
                const newCapacity = outputBufferCapacity * 2;
                gif_out_buffer_grow(outputBuffer, newCapacity, arena.pointer);
                outputBufferCapacity = newCapacity;
            }
        }

        gif_encoder_finish_frame(gifEncoder, outputBuffer);
        if (gif_out_buffer_capacity_left(outputBuffer) < 1024) {
            const newCapacity = outputBufferCapacity * 2;
            gif_out_buffer_grow(outputBuffer, newCapacity, arena.pointer);
            outputBufferCapacity = newCapacity;
        }
        gif_encoder_finish(gifEncoder, outputBuffer);

        const gifImage = arena.memory.getInt32(outputBuffer, true);
        const gifImageSize = arena.memory.getInt32(outputBuffer + 8, true);
        const gifImageArray = new Uint8Array(arena.memory.buffer, gifImage, gifImageSize);

        const gifBlob = new Blob([gifImageArray], { type: 'image/gif' });
        const gifUrl = URL.createObjectURL(gifBlob);

        arena.rewind(arenaSnapshot);

        return gifUrl;
    }

    return { encode };
})();

const canvas = new OffscreenCanvas(1, 1);
const context = canvas.getContext('2d', { willReadFrequently: true });

async function urlToImage(url) {
    const imageResponse = await fetch(url);
    const imageBlob = await imageResponse.blob();
    const image = await createImageBitmap(imageBlob);

    canvas.width = image.width;
    canvas.height = image.height;
    context.drawImage(image, 0, 0);

    // The order of bytes in ImageData is RGBA regardless of the platform endianness.
    const imageData = context.getImageData(0, 0, image.width, image.height);
    const components = Math.round(imageData.data.length / imageData.width / imageData.height);

    return {
        width: imageData.width,
        height: imageData.height,
        components,
        data: imageData.data,
    };
}

onmessage = async (event) => {
    const { type, payload } = event.data;

    if (type == 'encode') {
        const { imageUrl, config } = payload;

        const image = await urlToImage(imageUrl);
        if (image == null) {
            throw new Error(`Invalid image URL: ${imageUrl}`);
        }

        const gifEncoder = await gifEncoderPromise;

        const gifUrl = gifEncoder.encode(image, config);
        postMessage({
            type: 'result',
            payload: {
                gifUrl,
                width: image.width,
                height: image.height,
            },
        })
    }
};
