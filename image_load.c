#include "image_load.h"

#define STBI_NO_LINEAR
#define STBI_NO_HDR
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

bool image_load(char const *file_name, Image *image) {
    int image_width;
    int image_height;
    int image_bytes_per_pixel;
    u8 *image_data = stbi_load(file_name, &image_width, &image_height, &image_bytes_per_pixel, 0);
    if (image == NULL) {
        return false;
    }

    image->data = image_data;
    image->width = image_width;
    image->height = image_height;
    image->bytes_per_pixel = image_bytes_per_pixel;

    return true;
}

isize image_size(Image const *image) {
    return image->width * image->height * image->bytes_per_pixel;
}

void image_destroy(Image *image) {
    stbi_image_free(image->data);
}
