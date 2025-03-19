#ifndef IMAGE_LOAD_H
#define IMAGE_LOAD_H

#include <stdbool.h>

#include "common.h"

typedef struct {
    u8 *data;
    isize width;
    isize height;
    isize bytes_per_pixel;
} Image;

bool image_load(char const *file_name, Image *image);
isize image_size(Image const *image);
void image_destroy(Image *image);

#endif // IMAGE_LOAD_H
