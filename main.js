const ALIGNMENT = 16;
const BYTES_PER_PAGE = 64 * 1024;

class Arena {
    constructor() {
        this.currentPageCount = 16;
        this.maxPageCount = 8192;

        this.wasmMemory = new WebAssembly.Memory({
            initial: this.currentPageCount,
            maximum: this.maxPageCount,
        });
        this.memory = new DataView(this.wasmMemory.buffer);

        this.beginPointer = 8;
        this.memory.setInt32(this.beginPointer, 16, true);

        this.endPointer = 12;
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

        const pagesLeft = this.maxPageCount - this.currentPageCount;
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

        const padding = ~(this.begin + 1) & (ALIGNMENT - 1);
        this.ensure(padding + size);

        const newPointer = this.begin + padding;
        this.memory.setInt32(this.beginPointer, this.begin + padding + size, true);
        return newPointer;
    }
}

const arena = new Arena();

const wasmModule = await WebAssembly.instantiateStreaming(fetch('gif_encoder.wasm'), {
    env: {
        memory: arena.wasmMemory,
        reportAlloc: (arenaPointer, size) => {
            const padding = ~(arena.begin + 1) & (ALIGNMENT - 1);
            arena.ensure(padding + size);
        },

        pow: (base, exponent) => Math.pow(base, exponent),
        cbrt: (value) => Math.cbrt(value),
        round: (value) => Math.round(value),
    },
});

const {
    srgb_palette_black_and_white,
    srgb_palette_monochrome,
    srgb_palette_web_safe,

    srgb_to_float, srgb_to_linear, srgb_to_lab, srgb_to_oklab,
    float_to_srgb, linear_to_srgb, lab_to_srgb, oklab_to_srgb,

    colors_unique, colors_unique_inplace,

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
    gif_out_buffer_reset,

    gif_encoder_create,
    gif_encoder_start,
    gif_frame_delay,
    gif_encoder_finish,

    gif_encoder_start_frame,
    gif_encoder_feed_frame,
    gif_encode_whole_frame,
    gif_encoder_finish_frame,
} = wasmModule.instance.exports;

const canvas = new OffscreenCanvas(1, 1);
const context = canvas.getContext('2d', { willReadFrequently: true });
const imageInputElement = document.getElementById('image-input');

imageInputElement.addEventListener('input', async function() {
    const file = this.files[0];
    const fileUrl = URL.createObjectURL(file);

    const imageElement = await new Promise((resolve, reject) => {
        const imageElement = new Image();
        imageElement.onload = function() {
            resolve(this);
        };
        imageElement.onerror = function() {
            resolve(null);
        };
        imageElement.src = fileUrl;
    });
    if (imageElement == null) {
        return;
    }

    canvas.width = imageElement.width;
    canvas.height = imageElement.height;
    context.drawImage(imageElement, 0, 0);

    // The order of bytes in the ImageData is RGBA regardless of the platform endianness.
    const imageData = context.getImageData(0, 0, imageElement.width, imageElement.height);
    const components = Math.round(imageData.data.length / imageData.width / imageData.height);

    const image = {
        width: imageData.width,
        height: imageData.height,
        components,
        data: imageData.data,
    };
    processImage(image);
});

function processImage(image) {
    const arenaSnapshot = arena.snapshot();

    const components = image.components;

    const srgbImage = arena.alloc(image.width * image.height * components);
    {
        // Don't keep typed arrays in your hands for too long: as the memory grows they rot away.
        const srgbImageArray = new Uint8Array(
            arena.memory.buffer,
            srgbImage, image.width * image.height * components
        );
        srgbImageArray.set(image.data, 0);
    }

    const floatImage = srgb_to_float(
        srgbImage, image.width * image.height, components,
        arena.pointer
    );

    const colorCountPointer = arena.alloc(4);
    const srgbColors = srgb_palette_web_safe(
        colorCountPointer, components, arena.pointer
    );
    const colorCount = arena.memory.getInt32(colorCountPointer, true);
    const floatColors = srgb_to_float(
        srgbColors, colorCount, components,
        arena.pointer
    );

    const indexedImage = wasmModule.instance.exports.image_ordered_dither(
        floatImage, image.width, image.height, components,
        floatColors, colorCount,
        arena.pointer
    );

    let outputBufferCapacity = 64 * 1024;
    const outputBuffer = gif_out_buffer_create(outputBufferCapacity, arena.pointer);
    const gifEncoder = gif_encoder_create(arena.pointer);

    gif_encoder_start(
        gifEncoder,
        image.width, image.height,
        srgbColors, colorCount, components,
        outputBuffer
    );

    let indexedImageIter = indexedImage;
    const indexedImageEnd = indexedImage + image.width * image.height;

    gif_encoder_start_frame(gifEncoder, 0, 0, components, outputBuffer);

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

    const imageElement = new Image();
    imageElement.src = gifUrl;
    document.body.append(imageElement);

    arena.rewind(arenaSnapshot);
}
