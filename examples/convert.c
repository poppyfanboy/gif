// An example of converting an image into the GIF format using the library.
//
// convert -i <input image> -o <output GIF>
//        [-memory <arena size in MiB>]
//        [-color-space (srgb|linear|lab|oklab)]
//        [-palette (web-safe|monochrome|black-white|<file>)]
//        [-color-count <generated palette color count>]

#include <stdlib.h>     // malloc, free
#include <stddef.h>     // size_t, NULL
#include <stdio.h>      // FILE, fopen, fwrite, fclose, fprintf, stderr
#include <string.h>     // strcmp
#include <stdlib.h>     // strtoll
#include <errno.h>      // errno, ERANGE
#include <stdbool.h>    // bool, true, false

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "../src/gif_encoder.h"

typedef struct {
    u8 *begin;
    u8 *end;
} Arena;

typedef enum { SRGB, LINEAR, LAB, OKLAB } ColorSpace;

f32 *from_srgb(u8 const *srgb_colors, isize color_count, ColorSpace color_space, Arena *arena) {
    switch (color_space) {
    case SRGB:      return srgb_to_float(srgb_colors, color_count, arena);
    case LINEAR:    return srgb_to_linear(srgb_colors, color_count, arena);
    case LAB:       return srgb_to_lab(srgb_colors, color_count, arena);
    case OKLAB:     return srgb_to_oklab(srgb_colors, color_count, arena);
    }
}

u8 *into_srgb(f32 const *colors, isize color_count, ColorSpace color_space, Arena *arena) {
    switch (color_space) {
    case SRGB:      return float_to_srgb(colors, color_count, arena);
    case LINEAR:    return linear_to_srgb(colors, color_count, arena);
    case LAB:       return lab_to_srgb(colors, color_count, arena);
    case OKLAB:     return oklab_to_srgb(colors, color_count, arena);
    }
}

bool parse_int(char const *string, isize *integer) {
    char *string_end;
    *integer = (isize)strtoll(string, &string_end, 10);
    if (*string_end != '\0' || errno == ERANGE) {
        return false;
    }

    return true;
}

int main(int arg_count, char **args) {
    char const *input_file_name = NULL;
    char const *output_file_name = NULL;
    isize arena_capacity = 32 * 1024 * 1024;

    ColorSpace color_space = SRGB;
    enum { WEB_SAFE, MONOCHROME, BLACK_WHITE, GENERATE, FROM_FILE } palette = GENERATE;
    char const *palette_file_name = NULL;
    isize max_color_count = 256;

    for (isize i = 1; i < arg_count; i += 1) {
        if (strcmp(args[i], "-i") == 0 && i + 1 < arg_count) {
            input_file_name = args[++i];
        } else if (strcmp(args[i], "-o") == 0 && i + 1 < arg_count) {
            output_file_name = args[++i];
        } else if (strcmp(args[i], "-color-space") == 0 && i + 1 < arg_count) {
            char *color_space_arg = args[++i];

            if (strcmp(color_space_arg, "srgb") == 0) {
                color_space = SRGB;
            } else if (strcmp(color_space_arg, "linear") == 0) {
                color_space = LINEAR;
            } else if (strcmp(color_space_arg, "lab") == 0) {
                color_space = LAB;
            } else if (strcmp(color_space_arg, "oklab") == 0) {
                color_space = OKLAB;
            } else {
                fprintf(stderr, "Invalid color space: '%s'\n", color_space_arg);
                return 1;
            }
        } else if (strcmp(args[i], "-memory") == 0 && i + 1 < arg_count) {
            char *arena_capacity_arg = args[++i];

            if (!parse_int(arena_capacity_arg, &arena_capacity)) {
                fprintf(stderr, "Invalid arena capacity: '%s'\n", arena_capacity_arg);
                return 1;
            }
            arena_capacity *= 1024 * 1024;
        } else if (strcmp(args[i], "-palette") == 0 && i + 1 < arg_count) {
            char *palette_arg = args[++i];

            if (strcmp(palette_arg, "web-safe") == 0) {
                palette = WEB_SAFE;
            } else if (strcmp(palette_arg, "monochrome") == 0) {
                palette = MONOCHROME;
            } else if (strcmp(palette_arg, "black-white") == 0) {
                palette = BLACK_WHITE;
            } else {
                palette = FROM_FILE;
                palette_file_name = palette_arg;
            }
        } else if (strcmp(args[i], "-color-count") == 0) {
            char *color_count_arg = args[++i];

            if (
                !parse_int(color_count_arg, &max_color_count) ||
                !(1 <= max_color_count && max_color_count <= 256)
            ) {
                fprintf(stderr, "Invalid color count: '%s'\n", color_count_arg);
                return 1;
            }
        } else {
            fprintf(stderr, "Unknown or invalid option: '%s'\n", args[i]);
            return 1;
        }
    }
    if (input_file_name == NULL || output_file_name == NULL) {
        fprintf(stderr, "No input or output file names were provided.\n");
        return 1;
    }

    u8 *arena_memory = malloc(arena_capacity);
    if (arena_memory == NULL) {
        fprintf(stderr, "Failed to allocate memory for the arena.\n");
        return 1;
    }
    Arena arena = {arena_memory, arena_memory + arena_capacity};

    // Read the pixels from the file.

    f32 *pixels;
    int image_width, image_height;
    {
        u8 *srgb_pixels = stbi_load(input_file_name, &image_width, &image_height, NULL, 3);
        if (srgb_pixels == NULL) {
            fprintf(stderr, "Failed to load the input image: '%s'\n", input_file_name);
            return 1;
        }

        pixels = from_srgb(srgb_pixels, image_width * image_height, color_space, &arena);
        if (pixels == NULL) {
            fprintf(stderr, "Ran out of memory.\n");
            return 1;
        }
    }

    // Pick a color palette.

    f32 *colors;
    u8 *srgb_colors;
    isize color_count;
    if (palette == WEB_SAFE || palette == MONOCHROME || palette == BLACK_WHITE) {
        typedef u8 *(*PaletteGen)(isize *, void *);
        PaletteGen palette_gen = ((PaletteGen[]){
            [WEB_SAFE]      = srgb_palette_web_safe,
            [MONOCHROME]    = srgb_palette_monochrome,
            [BLACK_WHITE]   = srgb_palette_black_and_white,
        })[palette];

        srgb_colors = palette_gen(&color_count, &arena);
        colors = from_srgb(srgb_colors, color_count, color_space, &arena);
    } else if (palette == GENERATE) {
        colors = palette_by_median_cut(
            pixels, image_width * image_height,
            max_color_count,
            &color_count,
            &arena
        );
        srgb_colors = into_srgb(colors, color_count, color_space, &arena);
    } else if (palette == FROM_FILE) {
        int palette_width, palette_height;
        srgb_colors = stbi_load(palette_file_name, &palette_width, &palette_height, NULL, 3);
        if (srgb_colors == NULL) {
            fprintf(stderr, "Failed to load the file with the palette: '%s'\n", palette_file_name);
            return 1;
        }
        color_count = palette_width * palette_height;

        if (color_count == 0 || color_count > 256) {
            fprintf(stderr, "Too many colors in the palette: %td\n", color_count);
            return 1;
        }

        colors = from_srgb(srgb_colors, color_count, color_space, &arena);
    }
    if (colors == NULL || srgb_colors == NULL) {
        fprintf(stderr, "Ran out of memory.\n");
        return 1;
    }

    // Quantize the image using the picked color palette.

    GifColorIndex *indexed_pixels = image_quantize_for_gif(
        pixels, image_width * image_height,
        colors, color_count,
        &arena
    );
    if (indexed_pixels == NULL) {
        fprintf(stderr, "Ran out of memory.\n");
        return 1;
    }

    // Encode the GIF.

    FILE *output_file = fopen(output_file_name, "wb");
    if (output_file == NULL) {
        fprintf(stderr, "Failed to open an output file: '%s'\n", output_file_name);
        return 1;
    }

    GifEncoder *encoder = gif_encoder_create(&arena);
    GifOutputBuffer *out_buffer = gif_out_buffer_create(64 * 1024, &arena);
    if (encoder == NULL || out_buffer == NULL) {
        fprintf(stderr, "Ran out of memory.\n");
        return 1;
    }

    gif_encoder_start(encoder, image_width, image_height, srgb_colors, color_count, out_buffer);
    {
        GifColorIndex *indexed_pixel_iter = indexed_pixels;
        GifColorIndex *indexed_pixel_end = indexed_pixels + image_width * image_height;

        gif_encoder_start_frame(encoder, NULL, 0, out_buffer);

        while (indexed_pixel_iter != indexed_pixel_end) {
            indexed_pixel_iter += gif_encoder_feed_frame(
                encoder,
                indexed_pixel_iter,
                indexed_pixel_end - indexed_pixel_iter,
                out_buffer
            );
            if (gif_out_buffer_capacity_left(out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
                fwrite(out_buffer->data, 1, (size_t)out_buffer->encoded_size, output_file);
                gif_out_buffer_reset(out_buffer);
            }
        }

        gif_encoder_finish_frame(encoder, out_buffer);
        if (gif_out_buffer_capacity_left(out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
            fwrite(out_buffer->data, 1, (size_t)out_buffer->encoded_size, output_file);
            gif_out_buffer_reset(out_buffer);
        }
    }
    gif_encoder_finish(encoder, out_buffer);

    fwrite(out_buffer->data, 1, (size_t)out_buffer->encoded_size, output_file);
    fclose(output_file);

    return 0;
}
