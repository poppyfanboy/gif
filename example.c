#include <stdlib.h> // malloc, free
#include <stddef.h> // size_t, NULL
#include <stdio.h>  // FILE, fopen, fwrite, fclose

#include "common.h"
#include "image_load.h"
#include "gif_encoder.h"

#define ARENA_CAPACITY (64 * 1024 * 1024)
#define IO_BUFFER_SIZE (8 * 1024)

typedef GifArena Arena;
#define arena_alloc gif_arena_alloc

int main(void) {
    int exit_code = 0;

    u8 *arena_memory = NULL;
    FILE *input_file = NULL;
    FILE *output_file = NULL;

    arena_memory = malloc(ARENA_CAPACITY);
    if (arena_memory == NULL) {
        goto fail;
    }
    Arena arena = {arena_memory, arena_memory + ARENA_CAPACITY};

    // Load the input image into memory in its raw form.

    input_file = fopen("input.jpg", "rb");
    if (input_file == NULL) {
        goto fail;
    }

    isize image_width, image_height;
    if (!IMAGE_SIZE(input_file, &image_width, &image_height, arena)) {
        goto fail;
    }

    isize pixel_count = image_width * image_height;
    f32 *pixels = arena_alloc(&arena, pixel_count * 3 * sizeof(f32));
    {
        Arena arena_rewind = arena;

        u8 *pixels_srgb = IMAGE_LOAD(input_file, &arena);
        if (pixels_srgb == NULL) {
            goto fail;
        }
        srgb_to_linear(pixels_srgb, pixel_count, pixels);

        arena = arena_rewind;
    }

    // Pick a color palette and convert the input image, so that it only uses the selected colors.

    isize color_count = srgb_palette_web_safe(NULL);
    f32 *colors = arena_alloc(&arena, color_count * 3 * sizeof(f32));
    u8 *srgb_colors = arena_alloc(&arena, color_count * 3 * sizeof(u8));
    srgb_palette_web_safe(srgb_colors);
    srgb_to_linear(srgb_colors, color_count, colors);

    GifColorIndex *indexed_pixels = arena_alloc(&arena, pixel_count * sizeof(GifColorIndex));
    image_quantize_for_gif(pixels, pixel_count, colors, color_count, indexed_pixels, arena);

    // Encode the raw input image into the GIF format and write it into the output file.

    output_file = fopen("out.gif", "wb");
    if (output_file == NULL) {
        goto fail;
    };

    GifEncoder *encoder = gif_encoder_create(&arena);

    GifOutputBuffer out_buffer;
    gif_out_buffer_create(IO_BUFFER_SIZE, &out_buffer, &arena);

    gif_encoder_start(encoder, image_width, image_height, srgb_colors, color_count, &out_buffer);
    {
        GifColorIndex *indexed_pixel_iter = indexed_pixels;
        GifColorIndex *indexed_pixel_end = indexed_pixels + pixel_count;

        gif_encoder_start_frame(encoder, NULL, 0, &out_buffer);

        while (indexed_pixel_iter != indexed_pixel_end) {
            indexed_pixel_iter += gif_encoder_feed_frame(
                encoder,
                indexed_pixel_iter,
                indexed_pixel_end - indexed_pixel_iter,
                &out_buffer
            );
            if (gif_out_buffer_capacity_left(&out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
                fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
                gif_out_buffer_reset(&out_buffer);
            }
        }

        gif_encoder_finish_frame(encoder, &out_buffer);
        if (gif_out_buffer_capacity_left(&out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
            fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
            gif_out_buffer_reset(&out_buffer);
        }
    }
    gif_encoder_finish(encoder, &out_buffer);

    // Write any leftover buffered data into the file.
    fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);

    goto clean_up;

fail:
    exit_code = 1;
clean_up:
    if (output_file != NULL) {
        fclose(output_file);
    }
    if (input_file != NULL) {
        fclose(input_file);
    }
    free(arena_memory);

    return exit_code;
}
