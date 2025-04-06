#ifndef IMAGE_LOAD_H
#define IMAGE_LOAD_H

#include <stdio.h>      // FILE
#include <stdbool.h>    // bool, true, false

#include "common.h"

#define IMAGE_BYTES_PER_PIXEL 3

typedef struct {
    u8 *begin;
    u8 *end;
} ImageLoadArena;

isize memory_for_image_size_call(void);
bool image_size(FILE *file, isize *width, isize *height, ImageLoadArena scratch_arena);

isize max_memory_for_jpeg_image_load_call(isize width, isize height);
isize max_memory_for_png_image_load_call(FILE *file, isize width, isize height);

// For JPEG assumes an input image with 3 channels (8 bit) and 4:4:4 chroma subsampling.
// For PNG assumes an input image with 4 channels (8 bit), potentially with a palette.
isize max_memory_for_image_load_call(FILE *file, isize width, isize height);
u8 *image_load(FILE *file, ImageLoadArena *arena);

// Convenience wrapper macros that accept any arena-like type (in this case, one that has pointers
// to the beginning and the end of the available memory).

#define IMAGE_SIZE(file, width, height, scratch_arena) \
    image_size_impl(file, width, height, (scratch_arena).begin, (scratch_arena).end)

#define IMAGE_LOAD(file, arena) \
    image_load_impl(file, &(arena)->begin, &(arena)->end)

static inline bool image_size_impl(
    FILE *file,
    isize *width,
    isize *height,
    u8 *scratch_arena_begin,
    u8 *scratch_arena_end
) {
    ImageLoadArena scratch_arena = {.begin = scratch_arena_begin, .end = scratch_arena_end};
    return image_size(file, width, height, scratch_arena);
}

static inline u8 *image_load_impl(FILE *file, u8 **arena_begin, u8 **arena_end) {
    ImageLoadArena arena = {.begin = *arena_begin, .end = *arena_end};
    u8 *pixels = image_load(file, &arena);
    *arena_begin = arena.begin;
    *arena_end = arena.end;
    return pixels;
}

#endif // IMAGE_LOAD_H
