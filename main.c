// $ cc -DNDEBUG -O2 main.c
//
// GIF references:
// - https://www.w3.org/Graphics/GIF/spec-gif89a.txt
// - https://en.wikipedia.org/wiki/GIF

// TODO: Come up with a GIF encoder interface and extract it into a separate file.

#include <stdio.h>      // FILE, fopen, fclose, fwrite, fprintf, stderr
#include <stdint.h>     // uint8_t, uint16_t, ptrdiff_t
#include <stddef.h>     // size_t, NULL
#include <stdlib.h>     // malloc, free, exit, abort
#include <string.h>     // memmove, memcmp, memset
#include <stdbool.h>    // true, false, bool
#include <assert.h>     // assert
#include <math.h>       // powf, INFINITY

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define UNUSED(x) ((void) (x))

#define ARRAY_SIZE(array) ((isize) sizeof(array) / sizeof((array)[0]))

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef ptrdiff_t isize;
typedef float f32;

#define sizeof(x) ((isize) sizeof(x))

#define U8_WRITE(value, dest) ( \
    (dest)[0] = (u8) (value),   \
    (dest) + 1                  \
)

#define U16_WRITE_LE(value, dest) (             \
    (dest)[0] = (u8) ((u16) (value) & 0xff),    \
    (dest)[1] = (u8) ((u16) (value) >> 8),      \
    (dest) + 2                                  \
)

static inline u8 u8_min(u8 a, u8 b) {
    return a < b ? a : b;
}

static inline u8 u8_max(u8 a, u8 b) {
    return a > b ? a : b;
}

void file_write_all(FILE *file, u8 const *source, isize size) {
    isize bytes_written = (isize) fwrite(source, 1, (size_t) size, file);

    // > On success, fread() and fwrite() return the number of items read or written.
    // > This number equals the number of bytes transferred only when size is 1.
    // > If an error occurs, or the end of the file is reached, the return value is a short item
    // > count (or zero).
    //
    // Which means that (size != bytes_written) indicates an error while writing to the file.

    if (size != bytes_written) {
        // Failed to write to the output file.
        abort();
    }
}

// TODO: Use arenas for allocation?

void *alloc(isize size) {
    // > If size is 0, then malloc() returns a unique pointer value that can later be successfully
    // > passed to free().
    //
    // Return NULL expicitly instead, so that we would get a segfault, if we tried to dereference a
    // NULL pointer, as opposed to continuing business as usual, if the malloc happens to return a
    // valid pointer for a zero-sized allocation.

    if (size == 0) {
        return NULL;
    }

    void *memory = malloc((size_t) size);
    if (memory == NULL) {
        // Out of memory.
        abort();
    }

    return memory;
}

// Providing the size might help out an allocator, if it doesn't track allocation sizes itself (e.g.
// an arena allocator will be able to free an allocation in case it is located at the end of the
// buffer).

void dealloc(void *memory, isize size) {
    UNUSED(size);

    if (memory != NULL) {
        free(memory);
    }
}

#define GIF_MAX_WIDTH 0xffff
#define GIF_MAX_HEIGHT 0xffff

typedef u8 ColorIndex;

typedef struct {
    isize count;
    ColorIndex indices[];
} ColorSequence;

#define COLOR_ARRAY_MIN_CAPACITY 16

typedef struct {
    ColorIndex *indices;
    isize count;
    isize capacity;
} ColorArray;

void color_array_push(ColorArray *colors, ColorIndex index) {
    if (colors->count == colors->capacity) {
        isize new_capacity = colors->capacity * 2;
        if (new_capacity < COLOR_ARRAY_MIN_CAPACITY) {
            new_capacity = COLOR_ARRAY_MIN_CAPACITY;
        }

        ColorIndex *new_indices = alloc(new_capacity * sizeof(ColorIndex));
        memmove(new_indices, colors->indices, (size_t) (colors->count * sizeof(ColorIndex)));

        dealloc(colors->indices, colors->capacity * sizeof(ColorIndex));
        colors->indices = new_indices;
        colors->capacity = new_capacity;
    }

    colors->indices[colors->count] = index;
    colors->count += 1;
}

typedef isize LzwCode;
#define MAX_LZW_CODE 4095
#define LZW_CODE_MAX_BIT_LENGTH 12

typedef struct {
    isize count;
    ColorSequence *sequences[MAX_LZW_CODE + 1];
} LzwTable;

// TODO: Abstract away writing to a file. Write bytes to a dynamic array instead, and let the caller
// decide if they want to periodically flush the buffer to the file or keep stretching the buffer
// until it contains the entire GIF?

// Used to write "data sub-blocks" filled with LZW codes to the file.
typedef struct {
    FILE *output_file;

    u8 block[256];
    u8 *block_iter;

    // Next available bit position within the byte pointed to by the block_iter.
    isize bit_pos;
} LzwCodesWriter;

void lzw_writer_reset(LzwCodesWriter *writer) {
    // Skip over the size byte (we will set it later in the lzw_writer_flush function).
    writer->block_iter = writer->block + 1;
    writer->bit_pos = 0;

    memset(writer->block, 0, (size_t) sizeof(writer->block));
}

void lzw_writer_create(FILE *output_file, LzwCodesWriter *writer) {
    writer->output_file = output_file;
    lzw_writer_reset(writer);
}

bool lzw_writer_empty(LzwCodesWriter *writer) {
    return writer->block_iter == writer->block + 1 && writer->bit_pos == 0;
}

void lzw_writer_flush(LzwCodesWriter *writer) {
    // "Data sub-block" size including the size byte.
    isize block_size = writer->block_iter - writer->block;
    if (writer->bit_pos != 0) {
        block_size += 1;
    }

    assert(block_size > 0 && block_size - 1 <= 0xff);
    writer->block[0] = (u8) (block_size - 1);

    file_write_all(writer->output_file, writer->block, block_size);
    lzw_writer_reset(writer);
}

void lzw_write_code(LzwCodesWriter *writer, LzwCode code, isize code_bit_length) {
    assert(code_bit_length <= LZW_CODE_MAX_BIT_LENGTH);
    assert(code <= MAX_LZW_CODE);

    if (writer->block_iter - writer->block == sizeof(writer->block)) {
        lzw_writer_flush(writer);
    }

    isize code_bits_left = code_bit_length;
    LzwCode code_remainder = code;

    *writer->block_iter |= (u8) ((code_remainder << writer->bit_pos) & 0xff);

    if (code_bits_left + writer->bit_pos < 8) {
        writer->bit_pos += code_bits_left;
        code_bits_left = 0;
    } else {
        isize bits_written = 8 - writer->bit_pos;

        writer->bit_pos = 0;
        code_bits_left -= bits_written;
        code_remainder >>= bits_written;

        writer->block_iter += 1;
    }

    while (code_bits_left > 0) {
        assert(writer->bit_pos == 0);

        if (writer->block_iter - writer->block == sizeof(writer->block)) {
            lzw_writer_flush(writer);
        }

        *writer->block_iter |= (u8) (code_remainder & 0xff);

        if (code_bits_left < 8) {
            writer->bit_pos = code_bits_left;
            code_bits_left = 0;
        } else {
            isize bits_written = 8;

            writer->bit_pos = 0;
            code_bits_left -= bits_written;
            code_remainder >>= bits_written;

            writer->block_iter += 1;
        }
    }
}

#define RED_INDEX 0
#define GREEN_INDEX 1
#define BLUE_INDEX 2

typedef u32 RGB;

#define RGB(red, green, blue) ((RGB) red << 16 | (RGB) green << 8 | (RGB) blue)

typedef struct {
    isize count;
    RGB colors[256];
} ColorTable;

void color_table_black_and_white(ColorTable *color_table) {
    color_table->count = 2;
    color_table->colors[0] = 0x000000;
    color_table->colors[1] = 0xffffff;
}

void color_table_web_safe(ColorTable *color_table) {
    isize count = 0;
    for (RGB red = 0x00; red <= 0xff; red += 0x33) {
        for (RGB green = 0x00; green <= 0xff; green += 0x33) {
            for (RGB blue = 0x00; blue <= 0xff; blue += 0x33) {
                assert(count < ARRAY_SIZE(color_table->colors));
                color_table->colors[count] = RGB(red, green, blue);

                count += 1;
            }
        }
    }

    color_table->count = count;
}

#define GAMMA 2.2F

ColorIndex color_table_get_index(ColorTable *color_table, RGB needle) {
    f32 needle_red = powf((f32) ((needle >> 16) & 0xff) / 255.0F, GAMMA);
    f32 needle_green = powf((f32) ((needle >> 8) & 0xff) / 255.0F, GAMMA);
    f32 needle_blue = powf((f32) (needle & 0xff) / 255.0F, GAMMA);

    ColorIndex best_match = 0;
    f32 best_match_distance = INFINITY;

    RGB *color_iter = color_table->colors;
    while (color_iter < color_table->colors + color_table->count) {
        f32 iter_red = powf((f32) ((*color_iter >> 16) & 0xff) / 255.0F, GAMMA);
        f32 iter_green = powf((f32) ((*color_iter >> 8) & 0xff) / 255.0F, GAMMA);
        f32 iter_blue = powf((f32) (*color_iter & 0xff) / 255.0F, GAMMA);

        f32 iter_distance =
            (iter_red - needle_red) * (iter_red - needle_red) +
            (iter_green - needle_green) * (iter_green - needle_green) +
            (iter_blue - needle_blue) * (iter_blue - needle_blue);

        if (iter_distance < best_match_distance) {
            best_match = (ColorIndex) (color_iter - color_table->colors);
            best_match_distance = iter_distance;
        }

        color_iter += 1;
    }

    return best_match;
}

typedef struct {
    bool has_value;
    RGB color;
} ColorSetBucket;

typedef struct {
    ColorSetBucket *buckets;
    isize buckets_count;
    isize colors_count;
} ColorSet;

u32 rgb_hash(RGB value) {
    // Based on djb2 hashing function:
    // http://www.cse.yorku.ca/~oz/hash.html

    u32 hash = 5381;

    hash = ((hash << 5) + hash) + ((value >> 16) & 0xff);
    hash = ((hash << 5) + hash) + ((value >> 8) & 0xff);
    hash = ((hash << 5) + hash) + (value & 0xff);

    return hash;
}

void color_set_insert(ColorSet *colors, RGB new_color) {
    f32 fill_factor;
    if (colors->buckets_count == 0) {
        fill_factor = 1.0F;
    } else {
        fill_factor = (f32) colors->colors_count / colors->buckets_count;
    }

    if (fill_factor > 0.75F) {
        isize new_buckets_count = colors->buckets_count * 2;
        if (new_buckets_count < 16) {
            new_buckets_count = 16;
        }

        ColorSetBucket *new_buckets = malloc((size_t) (new_buckets_count * sizeof(ColorSetBucket)));
        memset(new_buckets, 0, (size_t) (new_buckets_count * sizeof(ColorSetBucket)));

        ColorSet rehashed_colors = {
            .buckets = new_buckets,
            .buckets_count = new_buckets_count,
            .colors_count = 0,
        };

        for (isize i = 0; i < colors->buckets_count; i += 1) {
            if (colors->buckets[i].has_value) {
                color_set_insert(&rehashed_colors, colors->buckets[i].color);
            }
        }

        free(colors->buckets);
        colors->buckets = rehashed_colors.buckets;
        colors->buckets_count = rehashed_colors.buckets_count;
        colors->colors_count = rehashed_colors.colors_count;
    }

    ColorSetBucket *bucket = &colors->buckets[rgb_hash(new_color) % colors->buckets_count];
    while (bucket->has_value) {
        if (bucket->color == new_color) {
            return;
        }

        bucket += 1;
        if (bucket == colors->buckets + colors->buckets_count) {
            bucket = colors->buckets;
        }
    }

    bucket->color = new_color;
    bucket->has_value = true;
    colors->colors_count += 1;
}

typedef struct {
    u8 *data;
    isize size;
    isize bytes_per_pixel;
} Image;

int rgb_compare_red_only(void const *left, void const *right) {
    u8 left_red = (*((RGB *) left) >> 16) & 0xff;
    u8 right_red = (*((RGB *) right) >> 16) & 0xff;

    return (int) left_red - (int) right_red;
}

int rgb_compare_green_only(void const *left, void const *right) {
    u8 left_green = (*((RGB *) left) >> 8) & 0xff;
    u8 right_green = (*((RGB *) right) >> 8) & 0xff;

    return (int) left_green - (int) right_green;
}

int rgb_compare_blue_only(void const *left, void const *right) {
    u8 left_blue = *((RGB *) left) & 0xff;
    u8 right_blue = *((RGB *) right) & 0xff;

    return (int) left_blue - (int) right_blue;
}

void median_cut_color_table(
    Image *image,
    isize target_colors_count,
    ColorTable *color_table
) {
    // FIXME: Add an assert that target_colors_count is less than 256

    ColorSet color_set = {NULL};

    u8 *image_data_iter = image->data;
    while (image_data_iter < image->data + image->size) {
        RGB color = RGB(
            image_data_iter[RED_INDEX],
            image_data_iter[GREEN_INDEX],
            image_data_iter[BLUE_INDEX]
        );
        color_set_insert(&color_set, color);

        image_data_iter += image->bytes_per_pixel;
    }

    RGB *colors = alloc(color_set.colors_count * sizeof(RGB));
    RGB *colors_end = colors + color_set.colors_count;

    {
        RGB *colors_iter = colors;
        ColorSetBucket *buckets_iter = color_set.buckets;
        while (buckets_iter < color_set.buckets + color_set.buckets_count) {
            if (buckets_iter->has_value) {
                assert(colors_iter != colors_end);

                *colors_iter = buckets_iter->color;
                colors_iter += 1;
            }

            buckets_iter += 1;
        }

        assert(colors_iter == colors_end);
    }

    isize segments_capacity = target_colors_count;

    // FIXME: Compute this from segment_begin_index and segment_end_index?
    isize segments_count = 1;

    RGB **segments_begins = alloc(segments_capacity * sizeof(RGB *));
    segments_begins[0] = colors;
    RGB **segments_ends = alloc(segments_capacity * sizeof(RGB *));
    segments_ends[0] = colors_end;

    isize segment_begin_index = 0;
    isize segment_end_index = 1;

    // FIXME: What if the number of target_colors_count is smaller than the actual nubmer of colors
    // in the image?

    while (segments_count < target_colors_count) {
        assert(segment_begin_index != segment_end_index);

        RGB *current_segment_begin = segments_begins[segment_begin_index];
        RGB *current_segment_end = segments_ends[segment_begin_index];

        // The image probably has less colors then we want in our palette, break out because
        // segments are naturally sorted by decreasing size, so the rest of the segments are goning
        // to be of size 1 as well.
        if (current_segment_end - current_segment_begin == 1) {
            break;
        }

        segment_begin_index = (segment_begin_index + 1) % segments_capacity;
        segments_count -= 1;

        u8 min_red = 0xff, min_green = 0xff, min_blue = 0xff;
        u8 max_red = 0x00, max_green = 0x00, max_blue = 0x00;
        {
            RGB *colors_iter = current_segment_begin;
            while (colors_iter != current_segment_end) {
                u8 red = (*colors_iter >> 16) & 0xff;
                u8 green = (*colors_iter >> 8) & 0xff;
                u8 blue = *colors_iter & 0xff;

                min_red = u8_min(min_red, red);
                min_green = u8_min(min_green, green);
                min_blue = u8_min(min_blue, blue);

                max_red = u8_max(max_red, red);
                max_green = u8_max(max_green, green);
                max_blue = u8_max(max_blue, blue);

                colors_iter += 1;
            }
        }

        int(*comparator)(void const *, void const *);

        if (max_red - min_red >= max_green - min_green && max_red - min_red >= max_blue - min_blue) {
            comparator = rgb_compare_red_only;
        } else if (max_green - min_green >= max_red - min_red && max_green - min_green >= max_blue - min_blue) {
            comparator = rgb_compare_green_only;
        } else {
            comparator = rgb_compare_blue_only;
        }

        assert(current_segment_begin < current_segment_end);
        qsort(
            current_segment_begin,
            current_segment_end - current_segment_begin,
            sizeof(RGB),
            comparator
        );

        isize current_segment_size = current_segment_end - current_segment_begin;
        assert(current_segment_size > 1);

        segments_begins[segment_end_index] = current_segment_begin;
        segments_ends[segment_end_index] = current_segment_begin + current_segment_size / 2;

        segments_begins[(segment_end_index + 1) % segments_capacity] = current_segment_begin + current_segment_size / 2;
        segments_ends[(segment_end_index + 1) % segments_capacity] = current_segment_end;

        segment_end_index = (segment_end_index + 2) % segments_capacity;
        segments_count += 2;
    }

    color_table->count = segments_count;
    RGB *palette_iter = color_table->colors;
    for (isize i = 0; i < segments_count; i += 1) {
        isize segment_index = (segment_begin_index + i) % segments_capacity;

        u32 red_sum = 0;
        u32 green_sum = 0;
        u32 blue_sum = 0;

        RGB *segment_iter = segments_begins[segment_index];
        RGB *segment_end = segments_ends[segment_index];
        while (segment_iter != segment_end) {
            // FIXME: Take gamma into account, maybe??

            red_sum += (*segment_iter >> 16) & 0xff;
            green_sum += (*segment_iter >> 8) & 0xff;
            blue_sum += *segment_iter & 0xff;

            segment_iter += 1;
        }

        isize colors_count = segments_ends[segment_index] - segments_begins[segment_index];
        RGB average = RGB(
            red_sum / colors_count,
            green_sum / colors_count,
            blue_sum / colors_count
        );

        *palette_iter = average;
        palette_iter += 1;
    }
}

int main(void) {
    int image_width;
    int image_height;
    int image_bytes_per_pixel;
    u8 *image = stbi_load("input.jpg", &image_width, &image_height, &image_bytes_per_pixel, 0);

    if (image == NULL) {
        // Failed to load an image.
        abort();
    }

    ColorTable color_table;

    median_cut_color_table(
        &(Image) {
            .data = image,
            .size = image_width * image_height * image_bytes_per_pixel,
            .bytes_per_pixel = image_bytes_per_pixel,
        },
        256,
        &color_table
    );

    if (image_width > GIF_MAX_WIDTH || image_height > GIF_MAX_HEIGHT) {
        // Image dimensions are too big for a GIF.
        abort();
    }

    FILE *output_file = fopen("out.gif", "wb");
    if (output_file == NULL) {
        // Failed to open an output file.
        abort();
    }


    // === GIF header ===

    // 6 bytes          for the "GIF89a" header
    // 7 bytes          for a logical screen descriptor
    // 256 * 3 bytes    for a max size global color table
    u8 header[6 + 7 + 256 * 3];

    u8 *header_iter = header;

    // GIF89a
    header_iter = U8_WRITE('G', header_iter);
    header_iter = U8_WRITE('I', header_iter);
    header_iter = U8_WRITE('F', header_iter);
    header_iter = U8_WRITE('8', header_iter);
    header_iter = U8_WRITE('9', header_iter);
    header_iter = U8_WRITE('a', header_iter);

    // Logical screen descriptor
    header_iter = U16_WRITE_LE(image_width, header_iter);
    header_iter = U16_WRITE_LE(image_height, header_iter);

    // FIXME: next_power_of_two is the worst variable name ever.
    assert(color_table.count > 0);
    int next_power_of_two =
        sizeof(unsigned int) * 8 -
        __builtin_clz((unsigned int) color_table.count);
    if ((color_table.count & (color_table.count - 1)) == 0) {
        next_power_of_two -= 1;
    }
    if (next_power_of_two < 1) {
        next_power_of_two = 1;
    }
    if (next_power_of_two > 8) {
        // color table is too big
        abort();
    }

    header_iter = U8_WRITE(0x80 | (u8) (next_power_of_two - 1), header_iter);

    header_iter = U8_WRITE(0, header_iter);
    header_iter = U8_WRITE(0, header_iter);

    // Global color table
    {
        RGB *color_iter = color_table.colors;

        for (isize i = 0; i < (1 << next_power_of_two); i += 1) {
            if (color_iter < color_table.colors + color_table.count) {
                header_iter = U8_WRITE(*color_iter >> 16, header_iter);
                header_iter = U8_WRITE((*color_iter >> 8) & 0xff, header_iter);
                header_iter = U8_WRITE(*color_iter & 0xff, header_iter);

                color_iter += 1;
            } else {
                header_iter = U8_WRITE(0x00, header_iter);
                header_iter = U8_WRITE(0x00, header_iter);
                header_iter = U8_WRITE(0x00, header_iter);
            }
        }
    }

    file_write_all(output_file, header, header_iter - header);


    // === Image header ===

    // 8 bytes  for a graphic control extension
    // 10 bytes for an image descriptor
    // 1 byte   for a minimum LZW code size
    u8 image_header[8 + 10 + 1];

    u8 *image_header_iter = image_header;

    // Graphic control extension
    image_header_iter = U8_WRITE(0x21, image_header_iter);
    image_header_iter = U8_WRITE(0xf9, image_header_iter);
    image_header_iter = U8_WRITE(4, image_header_iter);
    image_header_iter = U8_WRITE(0x08, image_header_iter);
    image_header_iter = U16_WRITE_LE(100, image_header_iter);
    image_header_iter = U8_WRITE(0, image_header_iter);
    image_header_iter = U8_WRITE(0, image_header_iter);

    // Image descriptor
    image_header_iter = U8_WRITE(0x2c, image_header_iter);
    image_header_iter = U16_WRITE_LE(0, image_header_iter);
    image_header_iter = U16_WRITE_LE(0, image_header_iter);
    image_header_iter = U16_WRITE_LE(image_width, image_header_iter);
    image_header_iter = U16_WRITE_LE(image_height, image_header_iter);
    image_header_iter = U8_WRITE(0x00, image_header_iter);

    // Minimum LZW code bit length
    // TODO: This value should be computed from the size of the color table.
    isize lzw_code_min_bit_length = next_power_of_two;
    *(image_header_iter++) = (u8) lzw_code_min_bit_length;

    file_write_all(output_file, image_header, image_header_iter - image_header);


    // === Image data ===

    LzwCodesWriter lzw_writer;
    lzw_writer_create(output_file, &lzw_writer);

    LzwCode lzw_clear_code = 1 << lzw_code_min_bit_length;
    LzwCode lzw_end_code = lzw_clear_code + 1;

    isize lzw_code_bit_length = lzw_code_min_bit_length + 1;
    LzwCode next_lzw_code = lzw_clear_code + 2;

    // Output a clear code because the spec says so:
    // > Encoders should output a Clear code as the first code of each image data stream.
    lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);

    LzwTable *lzw_table = alloc(sizeof(LzwTable));
    lzw_table->count = next_lzw_code;

    for (isize color_index = 0; color_index < (1 << next_power_of_two); color_index += 1) {
        ColorSequence *black_color_sequence = alloc(sizeof(ColorSequence) + sizeof(ColorIndex));
        black_color_sequence->count = 1;
        black_color_sequence->indices[0] = (ColorIndex) color_index;
        lzw_table->sequences[color_index] = black_color_sequence;
    }

    // Create a dummy sequence for the unused LZW codes, as well as for the "clear" and the "end of
    // information" codes, so that we don't have to handle these codes differently from others when
    // looking up sequences in the LZW table.
    ColorSequence *dummy_sequence = alloc(sizeof(ColorSequence));
    dummy_sequence->count = 0;

    for (
        LzwCode dummy_code = (1 << next_power_of_two);
        dummy_code < next_lzw_code;
        dummy_code += 1
    ) {
        lzw_table->sequences[dummy_code] = dummy_sequence;
    }

    u8 *image_iter = image;
    u8 *image_end = image + image_width * image_height * image_bytes_per_pixel;

    ColorArray current_sequence = {NULL};

    assert(image_end - image_iter > 0);

    // TODO: Abstract away looking up color indices from the color table.

    u8 red = image_iter[RED_INDEX];
    u8 green = image_iter[GREEN_INDEX];
    u8 blue = image_iter[BLUE_INDEX];
    image_iter += image_bytes_per_pixel;

    color_array_push(&current_sequence, color_table_get_index(&color_table, RGB(red, green, blue)));

    while (image_iter < image_end || current_sequence.count > 0) {
        bool nonexistent_sequence_found = false;
        LzwCode longest_sequence_lzw_code = 0;

        u8 *image_iter_rewind = image_iter;

        while (true) {
            bool sequence_already_exists = false;

            // TODO: Replace linear search with a hash table?

            for (LzwCode lzw_code = 0; lzw_code < lzw_table->count; lzw_code += 1) {
                if (current_sequence.count != lzw_table->sequences[lzw_code]->count) {
                    continue;
                }

                ColorSequence *sequence = lzw_table->sequences[lzw_code];

                // current_sequence is never empty at this point, so memcmp is not going to return 0
                // for dummy zero-length sequences.
                assert(current_sequence.count > 0);

                // BTW, apparently it is technically a UB to pass NULL to memcmp, as well as for
                // other string functions, even if the size is 0?
                //
                // https://stackoverflow.com/a/16363034
                // > Unless explicitly stated otherwise in the description of a particular function
                // > in this subclause, pointer arguments on such a call shall still have valid
                // > values

                if (memcmp(
                        current_sequence.indices,
                        sequence->indices,
                        (size_t) current_sequence.count
                    ) == 0
                ) {
                    sequence_already_exists = true;
                    longest_sequence_lzw_code = lzw_code;
                    break;
                }
            }

            if (!sequence_already_exists) {
                nonexistent_sequence_found = true;
                break;
            }

            if (image_iter < image_end) {
                u8 red = image_iter[RED_INDEX];
                u8 green = image_iter[GREEN_INDEX];
                u8 blue = image_iter[BLUE_INDEX];
                image_iter += image_bytes_per_pixel;

                color_array_push(
                    &current_sequence,
                    color_table_get_index(&color_table, RGB(red, green, blue))
                );
            } else {
                break;
            }
        }

        // > When the table is full, the encoder can chose to use the table as is, making no changes
        // > to it until the encoder chooses to clear it. The encoder during this time sends out
        // > codes that are of the maximum Code Size.
        //
        // Which means that you don't have to necessarily emit the clear code here?
        //
        // TODO: Simplify this part given the info above.

        if (nonexistent_sequence_found && lzw_table->count == MAX_LZW_CODE + 1) {
            lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);

            image_iter = image_iter_rewind;
            current_sequence.count = 1;

            for (LzwCode i = lzw_clear_code + 2; i < lzw_table->count; i += 1) {
                ColorSequence *sequence = lzw_table->sequences[i];
                dealloc(sequence, sizeof(ColorSequence) + sequence->count * sizeof(ColorIndex));
            }
            lzw_table->count = lzw_clear_code + 2;
            lzw_code_bit_length = lzw_code_min_bit_length + 1;

            continue;
        }

        lzw_write_code(&lzw_writer, longest_sequence_lzw_code, lzw_code_bit_length);

        if (nonexistent_sequence_found) {
            // In here current_sequence length is at least 2 because all possible sequences of
            // length 1 are guaranteed to exist in the table.


            if (lzw_table->count != MAX_LZW_CODE + 1) {
                ColorSequence *new_sequence = alloc(
                    sizeof(ColorSequence) + current_sequence.count * sizeof(ColorIndex)
                );
                new_sequence->count = current_sequence.count;
                memmove(
                    new_sequence->indices,
                    current_sequence.indices,
                    (size_t) (current_sequence.count * sizeof(ColorIndex))
                );

                lzw_table->sequences[lzw_table->count] = new_sequence;

                // Increase the LZW code length here because the spec says so:
                // > Whenever the LZW code value would exceed the current code length,
                // > the code length is increased by one. The packing/unpacking of these
                // > codes must then be altered to reflect the new code length.

                if ((lzw_table->count & (lzw_table->count - 1)) == 0) {
                    lzw_code_bit_length += 1;
                }

                lzw_table->count += 1;
            }

            ColorIndex next_after_existing = current_sequence.indices[current_sequence.count - 1];
            current_sequence.count = 1;
            current_sequence.indices[0] = next_after_existing;
        } else {
            assert(image_iter == image_end);
            current_sequence.count = 0;
        }
    }

    lzw_write_code(&lzw_writer, lzw_end_code, lzw_code_bit_length);

    if (!lzw_writer_empty(&lzw_writer)) {
        lzw_writer_flush(&lzw_writer);
    }

    // TODO: Make it possible to encode multiple images into a single GIF.

    u8 block_terminator[] = {0x00};
    file_write_all(output_file, block_terminator, sizeof(block_terminator));


    // === GIF trailer ===

    u8 trailer[] = {0x3b};
    file_write_all(output_file, trailer, sizeof(trailer));


    // Ensure that the FILE buffer is flushed.
    fclose(output_file);

    return 0;
}
