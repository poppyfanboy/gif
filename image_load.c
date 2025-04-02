#include "image_load.h"

#include <stddef.h> // size_t, NULL
#include <string.h> // memcpy
#include <assert.h> // assert

#define ARENA_ALIGNMENT 16

static ImageLoadArena *current_arena = NULL;

static void *arena_alloc(size_t size) {
    assert(current_arena != NULL);

    isize padding = (-((uptr) current_arena->begin)) & (ARENA_ALIGNMENT - 1);
    isize memory_left = current_arena->end - current_arena->begin - padding;
    if (memory_left < 0 || (size_t) memory_left < size) {
        // stbi_load will exit gracefully and return NULL, if we fail to allocate memory.
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
        // (Doesn't seem to happen in practice with stbi_load.)
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

u8 *image_load(char const *file_name, isize *width_out, isize *height_out, ImageLoadArena *arena) {
    current_arena = arena;
    ImageLoadArena arena_checkpoint = *arena;

    int width = 0;
    int height = 0;
    int channel_count;
    u8 *pixels = stbi_load(file_name, &width, &height, &channel_count, IMAGE_BYTES_PER_PIXEL);
    if (pixels == NULL) {
        *arena = arena_checkpoint;
    }

    *width_out = width;
    *height_out = height;

    current_arena = NULL;
    return pixels;
}
