#include "image_load.h"

#include <stddef.h> // size_t, NULL
#include <string.h> // memcpy
#include <assert.h> // assert
#include <stdio.h>  // fseek, ftell

#define ARENA_ALIGNMENT 16

static ImageLoadArena *current_arena = NULL;

static void *arena_alloc(size_t size) {
    assert(current_arena != NULL);

    isize padding = (-((uptr) current_arena->begin)) & (ARENA_ALIGNMENT - 1);
    isize memory_left = current_arena->end - current_arena->begin - padding;
    if (memory_left < 0 || (size_t) memory_left < size) {
        // stb_image will fail gracefully and return NULL or false, if we run out of memory.
        return NULL;
    }

    void *ptr = current_arena->begin + padding;
    current_arena->begin += padding + (isize) size;
    return ptr;
}

static void arena_free(void *ptr) {
    assert(current_arena != NULL);

    // Cannot do anything special here, because stb_image does not give us the allocation size.
    (void) ptr;
}

static void *arena_realloc(void *old_ptr, size_t old_size, size_t new_size) {
    assert(current_arena != NULL);

    if ((u8 *) old_ptr + old_size < current_arena->begin) {
        void *new_ptr = arena_alloc(new_size);
        memcpy(new_ptr, old_ptr, old_size);
        return new_ptr;
    } else {
        // Try to extend an old allocation, if it happens to be the topmost one.
        //
        // This actually saves up a little bit of memory when decoding a PNG, specifically while
        // concatenating all IDAT chunks into a single dynamic array.
        current_arena->begin -= old_size;
        return arena_alloc(new_size);
    }
}

#define STBI_MALLOC(size) arena_alloc(size)
#define STBI_REALLOC_SIZED(ptr, old_size, new_size) arena_realloc(ptr, old_size, new_size)
#define STBI_FREE(ptr) arena_free(ptr)

#define STBI_NO_LINEAR
#define STBI_ONLY_PNG
#define STBI_ONLY_JPEG
#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static inline isize isize_max(isize lhs, isize rhs) {
    return lhs > rhs ? lhs : rhs;
}

isize memory_for_image_size_call(void) {
    isize memory_for_jpeg = ARENA_ALIGNMENT + sizeof(stbi__jpeg);

    // JPEG is here as well, because stb_image first tries to interpret a file as a JPEG image.
    isize memory_for_png = ARENA_ALIGNMENT + sizeof(stbi__jpeg);

    return isize_max(memory_for_jpeg, memory_for_png);
}

bool image_size(FILE *file, isize *width_out, isize *height_out, ImageLoadArena scratch_arena) {
    current_arena = &scratch_arena;

    int width = 0;
    int height = 0;
    int channel_count = 0;
    bool success = stbi_info_from_file(file, &width, &height, &channel_count) == 1;

    *width_out = width;
    *height_out = height;

    current_arena = NULL;
    return success;
}

isize max_memory_for_jpeg_image_load_call(isize width, isize height) {
    isize memory = 0;

    // Brutforce over image formats, including JPEG.
    memory += ARENA_ALIGNMENT + sizeof(stbi__jpeg);

    // Decode the JPEG header (the struct allocated above is not getting reused).
    memory += ARENA_ALIGNMENT + sizeof(stbi__jpeg);

    isize channel_count = 3;

    // Worst case scenario is 4:4:4 chroma subsampling (1x1, 1x1, 1x1).
    // Use 4:2:0 (2x2, 1x1, 1x1) or 4:2:2 (2x1, 1x1, 1x1) instead here, as they are more common?
    isize subsample_x[]  = {1, 1, 1};
    isize subsample_y[] = {1, 1, 1};
    isize max_subsample_x =
        isize_max(isize_max(subsample_x[0], subsample_x[1]), subsample_x[2]);
    isize max_subsample_y =
        isize_max(isize_max(subsample_y[0], subsample_y[1]), subsample_y[2]);

    isize mcu_width = max_subsample_y * 8;
    isize mcu_count_x = (width + mcu_width - 1) / mcu_width;
    isize mcu_height = max_subsample_x * 8;
    isize mcu_count_y = (height + mcu_height - 1) / mcu_height;

    // Allocate this buffer for each channel (see the stbi__process_frame_header function).
    for (isize channel = 0; channel < channel_count; channel += 1) {
        memory += ARENA_ALIGNMENT + (
            mcu_count_x * subsample_x[channel] * 8 *
            mcu_count_y * subsample_y[channel] * 8
        ) + 15;
    }

    // A line buffer per channel (see the load_jpeg_image function).
    memory += (ARENA_ALIGNMENT + width + 3) * channel_count;

    // A buffer for the decoded image.
    memory += ARENA_ALIGNMENT + width * height * channel_count;

    return memory;
}

isize max_memory_for_png_image_load_call(FILE *file, isize width, isize height) {
    isize memory = 0;

    // Brutforce over image formats, yet again, including JPEG.
    memory += ARENA_ALIGNMENT + sizeof(stbi__jpeg);

    // stb_image concatenates data from all IDAT chunks into a single dynamic array.
    // Estimate the final size of the array to be the same as the size of the input image.
    //
    // In the worst case stb_image runs out of the array capacity right at the end of the last
    // IDAT chunk, as a result, we will allocate twice the amount of memory we actually needed
    // (the growth factor of the array is 2).
    long file_current_pos = ftell(file);
    fseek(file, 0, SEEK_END);
    long file_end_pos = ftell(file);
    fseek(file, file_current_pos, SEEK_SET);
    memory += ARENA_ALIGNMENT + 2 * (file_end_pos - file_current_pos);

    isize channel_count = 4;
    isize depth = 8;

    isize bytes_per_line = (width * depth + 7) / 8;
    // A buffer for the inflated ZLIB data (see the IEND case within the stbi__parse_png_file)
    // in case we are dealing with a non-paletted PNG image.
    //
    // An additional "+ height" is for the filter mode per row.
    memory += ARENA_ALIGNMENT + bytes_per_line * height * (channel_count + 1) + height;

    // This is the same ZLIB decoded buffer, but in case the input PNG uses an 8 bit palette.
    //
    // For paletted PNGs stb_image also allocates an additional buffer to expand color indices
    // into actual colors (see stbi__expand_png_palette), but we kind of already accounted for
    // this with the allocation above.
    memory += ARENA_ALIGNMENT + width * height + height;

    // A buffer for the raw decoded image.
    memory += ARENA_ALIGNMENT + width * height * (depth == 16 ? 2 : 1) * channel_count;

    // Two scanlines for filtering (see the stbi__create_png_image_raw function).
    memory += ARENA_ALIGNMENT + (width * depth * channel_count + 7) / 8 * 2;

    // Yet another buffer for the raw decoded image, but with the exact requested channel count
    // which is 3 (see the stbi__do_png function).
    if (channel_count != 3) {
        memory += ARENA_ALIGNMENT + width * height * channel_count;
    }

    return memory;
}

isize max_memory_for_image_load_call(FILE *file, isize width, isize height) {
    return isize_max(
        max_memory_for_jpeg_image_load_call(width, height),
        max_memory_for_png_image_load_call(file, width, height)
    );
}

u8 *image_load(FILE *file, ImageLoadArena *arena) {
    current_arena = arena;
    ImageLoadArena arena_checkpoint = *arena;

    int width;
    int height;
    int channel_count;
    u8 *pixels = stbi_load_from_file(file, &width, &height, &channel_count, IMAGE_BYTES_PER_PIXEL);
    if (pixels == NULL) {
        *arena = arena_checkpoint;
    }

    current_arena = NULL;
    return pixels;
}
