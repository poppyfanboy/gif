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

const wasmMemory = new WebAssembly.Memory({ initial: INITIAL_PAGE_COUNT, maximum: MAX_PAGE_COUNT });

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

async function fileToImage(file) {
    const fileUrl = URL.createObjectURL(file);

    const image = await new Promise((resolve) => {
        const image = new Image();
        image.onload = function() {
            resolve(this);
        };
        image.onerror = function() {
            resolve(null);
        };
        image.src = fileUrl;
    });

    if (image == null) {
        return null;
    }

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
        element: image,
        data: imageData.data,
    };
}

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
function imageToGif(image, config) {
    const arenaSnapshot = arena.snapshot();

    const srgbImage = arena.alloc(image.width * image.height * image.components);
    {
        // Don't keep typed arrays in your hands for too long: as the memory grows they rot away.
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

    const imageElement = new Image();
    imageElement.src = gifUrl;
    imageElement.width = image.width;
    imageElement.height = image.height;

    arena.rewind(arenaSnapshot);

    return imageElement;
}

document.addEventListener('dragover', (event) => {
    event.preventDefault();
});

document.addEventListener('drop', (event) => {
    event.preventDefault();
});

class ImageInput extends HTMLElement {
    constructor() {
        super();

        const shadowRoot = this.attachShadow({ mode: 'open' });
        const template = document.getElementById('image-input-template')
        shadowRoot.appendChild(template.content.cloneNode(true));
    }

    connectedCallback() {
        const overlay = this.shadowRoot.querySelector('.overlay');

        document.addEventListener('drop', () => {
            overlay.classList.remove('overlay--visible');
            overlay.classList.remove('overlay--active');
        });

        document.addEventListener('dragenter', (event) => {
            if (event.relatedTarget == null) {
                overlay.classList.add('overlay--visible');
            }
        });

        document.addEventListener('dragleave', (event) => {
            if (event.relatedTarget == null) {
                overlay.classList.remove('overlay--visible');
                overlay.classList.remove('overlay--active');
            }
        });

        const dropZone = this.shadowRoot.querySelector('.drop-zone');

        dropZone.addEventListener('dragenter', () => {
            overlay.classList.add('overlay--active');
        });

        dropZone.addEventListener('dragleave', (event) => {
            if (event.relatedTarget == null || !dropZone.contains(event.relatedTarget)) {
                overlay.classList.remove('overlay--active');
            }
        });

        const fileInput = this.shadowRoot.querySelector('input');

        dropZone.addEventListener('drop', async (event) => {
            const droppedItem = event.dataTransfer.items[0];
            if (droppedItem.kind != 'file') {
                return;
            }

            const image = await fileToImage(droppedItem.getAsFile());
            if (image == null) {
                return;
            }

            this.showUploadedImage(image.element);

            const changeEvent = new CustomEvent('change', {
                composed: true,
                detail: { image },
            });
            this.shadowRoot.dispatchEvent(changeEvent);
        });

        fileInput.addEventListener('change', async () => {
            if (fileInput.files.length == 0) {
                return;
            }

            const image = await fileToImage(fileInput.files[0]);
            if (image == null) {
                return;
            }

            fileInput.value = '';
            this.showUploadedImage(image.element);

            const changeEvent = new CustomEvent('change', {
                composed: true,
                detail: { image },
            });
            this.shadowRoot.dispatchEvent(changeEvent);
        });
    }

    showUploadedImage(imageElement) {
        const preview = this.shadowRoot.querySelector('.preview');
        preview.style.backgroundImage = `url('${imageElement.src}')`;

        const dropZone = this.shadowRoot.querySelector('.drop-zone');
        dropZone.style.aspectRatio = `${imageElement.width} / ${imageElement.height}`;
    }
}

customElements.define('image-input', ImageInput);

class GifConfig {
    /** @type {PaletteType} */
    paletteType = 'fixed';

    constructor() {
        this.colorSpaceSelect = document.getElementById('color-space-select');

        this.ditheringAlgorithmSelect = document.getElementById('dithering-algorithm-select');

        const fixedPaletteLabel = document.getElementById('fixed-palette-label');
        this.paletteSelect = document.getElementById('fixed-palette-select');

        const paletteAlgorthmLabel = document.getElementById('palette-algorithm-label');
        this.paletteAlgorithmSelect = document.getElementById('palette-algorithm-select');

        const paletteColorCountLabel = document.getElementById('palette-color-count-label');
        this.paletteColorCountInput = document.getElementById('palette-color-count-input');

        const fixedPaletteRadio = document.querySelector(
            '[name=palette-type][value=fixed]'
        );
        fixedPaletteRadio.checked = true;
        fixedPaletteRadio.addEventListener('change', (event) => {
            if (event.target.checked) {
                this.paletteType = 'fixed';

                fixedPaletteLabel.classList.remove('gif-config_label--hidden');

                paletteAlgorthmLabel.classList.add('gif-config_label--hidden');
                paletteColorCountLabel.classList.add('gif-config_label--hidden');
            }
        });

        const generatedPaletteRadio = document.querySelector(
            '[name=palette-type][value=generated]'
        );
        generatedPaletteRadio.addEventListener('change', (event) => {
            if (event.target.checked) {
                this.paletteType = 'generated';

                paletteAlgorthmLabel.classList.remove('gif-config_label--hidden');
                paletteColorCountLabel.classList.remove('gif-config_label--hidden');

                fixedPaletteLabel.classList.add('gif-config_label--hidden');
            }
        });
    }

    /** @returns {FixedPalette} */
    get fixedPalette() {
        return this.paletteSelect.value;
    }

    /** @returns {ColorSpace} */
    get colorSpace() {
        return this.colorSpaceSelect.value;
    }

    /** @returns {DitheringAlgorithm} */
    get ditheringAlgorithm() {
        return this.ditheringAlgorithmSelect.value;
    }

    /** @returns {PaletteAlgorithm} */
    get paletteAlgorithm() {
        return this.paletteAlgorithmSelect.value;
    }

    /** @returns {number} */
    get paletteColorCount() {
        const colorCount = Number.parseInt(this.paletteColorCountInput.value);

        if (Number.isNaN(colorCount)) {
            return 16;
        } else {
            return colorCount;
        }
    }
}

class App {
    inputImage = null;
    convertedImage = null;

    constructor() {
        const imageInput = document.querySelector('image-input');
        imageInput.addEventListener('change', (event) => {
            const { image } = event.detail;

            this.inputImage = image;
            this.convertedImage = null;

            const imageOutput = document.getElementById('image-output');
            imageOutput.style.aspectRatio = `${image.width} / ${image.height}`;

            const imageOutputPreview = document.getElementById('image-output-preview');
            imageOutputPreview.style.backgroundImage = null;
        });

        const gifConfig = new GifConfig();

        const convertButton = document.getElementById('convert-button');
        convertButton.addEventListener('click', async () => {
            if (this.inputImage == null) {
                return;
            }

            this.convertedImage = imageToGif(this.inputImage, gifConfig);
            this.showConvertedImage(this.convertedImage);
        });

        const downloadButton = document.getElementById('download-button');
        downloadButton.addEventListener('click', () => {
            if (this.convertedImage == null) {
                return;
            }

            const anchorElement = document.createElement('a');
            anchorElement.href = this.convertedImage.src;
            anchorElement.download = 'image.gif';
            anchorElement.click();
        });
    }

    showConvertedImage(image) {
        const preview = document.getElementById('image-output-preview');
        const imageOutput = document.getElementById('image-output');

        preview.style.backgroundImage = `url(${image.src})`;
        imageOutput.style.aspectRatio = `${image.width} / ${image.height}`;
    }
}

const app = new App();
