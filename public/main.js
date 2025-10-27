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
        URL.revokeObjectURL(fileUrl);
        return null;
    }

    return image;
}

document.addEventListener('dragover', (event) => {
    event.preventDefault();
});

document.addEventListener('drop', (event) => {
    event.preventDefault();
});

class ImageInput extends HTMLElement {
    static observedAttributes = ['disabled'];

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

            this.showUploadedImage(image);

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

            this.showUploadedImage(image);

            const changeEvent = new CustomEvent('change', {
                composed: true,
                detail: { image },
            });
            this.shadowRoot.dispatchEvent(changeEvent);
        });
    }

    attributeChangedCallback(name, _oldValue, newValue) {
        if (name == 'disabled') {
            const input = this.shadowRoot.querySelector('input');

            if (newValue == 'true') {
                input.disabled = true;
            } else if (newValue == 'false') {
                input.disabled = false;
            }
        }
    }

    showUploadedImage(image) {
        const preview = this.shadowRoot.querySelector('.preview');
        preview.style.backgroundImage = `url('${image.src}')`;

        const dropZone = this.shadowRoot.querySelector('.drop-zone');
        dropZone.style.aspectRatio = `${image.width} / ${image.height}`;
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

    toObject() {
        return {
            paletteType: this.paletteType,
            fixedPalette: this.fixedPalette,
            colorSpace: this.colorSpace,
            ditheringAlgorithm: this.ditheringAlgorithm,
            paletteAlgorithm: this.paletteAlgorithm,
            paletteColorCount: this.paletteColorCount,
        };
    }
}

class App {
    inputImageUrl = null;
    convertedImageUrl = null;

    constructor() {
        this.gifWorker = new Worker('gif-worker.js');
        this.gifConfig = new GifConfig();

        const imageInput = document.querySelector('image-input');
        imageInput.addEventListener('change', (event) => {
            const { image } = event.detail;

            if (this.inputImageUrl != null) {
                URL.revokeObjectURL(this.inputImageUrl);
            }
            this.inputImageUrl = image.src;

            if (this.convertedImageUrl != null) {
                URL.revokeObjectURL(this.convertedImageUrl);
            }

            const imageOutput = document.getElementById('image-output');
            imageOutput.style.aspectRatio = `${image.width} / ${image.height}`;

            const imageOutputPreview = document.getElementById('image-output-preview');
            imageOutputPreview.style.backgroundImage = null;
        });

        const convertButton = document.getElementById('convert-button');
        convertButton.addEventListener('click', async () => {
            if (this.inputImageUrl == null) {
                return;
            }

            this.gifWorker.postMessage({
                type: 'encode',
                payload: {
                    imageUrl: this.inputImageUrl,
                    config: this.gifConfig.toObject(),
                }
            });

            convertButton.disabled = true;
            convertButton.classList.add('loading');
            imageInput.setAttribute('disabled', true);
        });

        this.gifWorker.onmessage = (event) => {
            const { type, payload } = event.data;

            if (type == 'result') {
                const { gifUrl, width, height } = payload;

                const preview = document.getElementById('image-output-preview');
                const imageOutput = document.getElementById('image-output');

                preview.style.backgroundImage = `url(${gifUrl})`;
                imageOutput.style.aspectRatio = `${width} / ${height}`;

                this.convertedImageUrl = gifUrl;

                convertButton.disabled = false;
                convertButton.classList.remove('loading');
                imageInput.setAttribute('disabled', false);
            }
        };

        this.gifWorker.onerror = () => {
            convertButton.disabled = false;
            convertButton.classList.remove('loading');
            imageInput.setAttribute('disabled', false);
        };

        const downloadButton = document.getElementById('download-button');
        downloadButton.addEventListener('click', () => {
            if (this.convertedImageUrl == null) {
                return;
            }

            const anchorElement = document.createElement('a');
            anchorElement.href = this.convertedImageUrl;
            anchorElement.download = 'image.gif';
            anchorElement.click();
        });
    }
}

const app = new App();
