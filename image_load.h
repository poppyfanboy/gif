#ifndef IMAGE_LOAD_H
#define IMAGE_LOAD_H

#include "common.h"

#define IMAGE_BYTES_PER_PIXEL 3

typedef struct {
    u8 *begin;
    u8 *end;
} ImageLoadArena;

u8 *image_load(char const *file_name, isize *width, isize *height, ImageLoadArena *arena);

// A convenience wrapper macro that accepts a pointer to any arena-like type (in this case, one that
// has pointers to the beginning and the end of the available memory).
#define IMAGE_LOAD(file_name, width, height, arena) \
    image_load_impl(file_name, width, height, &(arena)->begin, &(arena)->end)

static inline u8 *image_load_impl(
    char const *file_name,
    isize *width,
    isize *height,
    u8 **arena_begin,
    u8 **arena_end
) {
    ImageLoadArena arena = {.begin = *arena_begin, .end = *arena_end};
    u8 *pixels = image_load(file_name, width, height, &arena);
    *arena_begin = arena.begin;
    *arena_end = arena.end;
    return pixels;
}

#endif // IMAGE_LOAD_H
