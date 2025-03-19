// $ cc -DNDEBUG -O2 -c image_load.c -o image_load.o
// $ cc -DNDEBUG -O2 -lm image_load.o main.c
//
// GIF references:
// - https://www.w3.org/Graphics/GIF/spec-gif89a.txt
// - https://en.wikipedia.org/wiki/GIF

// TODO: Come up with a GIF encoder interface and extract it into a separate file.

#include <stdio.h>      // FILE, fopen, fclose, fwrite, stderr, fprintf
#include <stddef.h>     // size_t, NULL
#include <stdlib.h>     // malloc, realloc, free, abort
#include <string.h>     // memmove, memcmp, memset
#include <stdbool.h>    // true, false, bool
#include <assert.h>     // assert
#include <math.h>       // powf, INFINITY
#include <stdarg.h>     // va_list, va_start, va_end, va_arg

#include "common.h"
#include "image_load.h"

static inline u8 *u8_write(u8 value, u8 *dest) {
    dest[0] = value;
    return dest + 1;
}

static inline u8 *u16_write_le(u16 value, u8 *dest) {
    dest[0] = (u8) (value & 0xff);
    dest[1] = (u8) (value >> 8);
    return dest + 2;
}

static inline u8 u8_min(u8 lhs, u8 rhs) {
    return lhs < rhs ? lhs : rhs;
}

static inline u8 u8_max(u8 lhs, u8 rhs) {
    return lhs > rhs ? lhs : rhs;
}

static inline i32 u32_clz(u32 value) {
#if defined(__clang__) || defined(__GNUC__)
    return __builtin_clz(value);
#elif defined(_MSC_VER)
    return (i32) __lzcnt(value);
#else
    i32 result = 0;
    if ((value & 0xffff0000) == 0) {
        result += 16;
        value <<= 16;
    }
    if ((value & 0xff000000) == 0) {
        result += 8;
        value <<= 8;
    }
    if ((value & 0xf0000000) == 0) {
        result += 4;
        value <<= 4;
    }
    if ((value & 0xc0000000) == 0) {
        result += 2;
        value <<= 2;
    }
    if ((value & 0x80000000) == 0) {
        result += 1;
    }
    return result;
#endif
}

i32 i32_log2_ceil(i32 value) {
    assert(value > 0);

    i32 log2_ceil = sizeof(u32) * 8 - u32_clz((u32) value);
    if ((value & (value - 1)) == 0) {
        log2_ceil -= 1;
    }

    return log2_ceil;
}

typedef enum {
    GIF_OUTPUT_IMAGE_WRITE_ERROR,
    GIF_OUT_OF_MEMORY,
    GIF_INPUT_IMAGE_READ_ERROR,
    GIF_ERROR_INPUT_IMAGE_TOO_BIG,
    GIF_ERROR_COLOR_TABLE_OVERFLOW,
} GifErrorKind;

// Any additional info associated with the error is passed through the variadic parameters.
void gif_error_callback_impl(GifErrorKind error_kind, ...) {
    va_list error_data;

    switch (error_kind) {
    case GIF_OUTPUT_IMAGE_WRITE_ERROR: {
        fprintf(stderr, "[ERROR] Failed to write the output image.\n");
    } break;

    case GIF_OUT_OF_MEMORY: {
        fprintf(stderr, "[ERROR] Out of memory.\n");
    } break;

    case GIF_INPUT_IMAGE_READ_ERROR: {
        fprintf(stderr, "[ERROR] Failed to read the input image.\n");
    } break;

    case GIF_ERROR_INPUT_IMAGE_TOO_BIG: {
        va_start(error_data, error_kind);
        isize input_image_width = va_arg(error_data, isize);
        isize input_image_height = va_arg(error_data, isize);

        fprintf(
            stderr,
            "[ERROR] Input image dimensions ("FMT_ISIZE"x"FMT_ISIZE") are too big for a GIF.\n",
            input_image_width,
            input_image_height
        );

        va_end(error_data);
    } break;

    case GIF_ERROR_COLOR_TABLE_OVERFLOW: {
        fprintf(stderr, "[ERROR] GIF color table overflow.\n");
    } break;
    }
}

#ifdef NDEBUG
    #define GIF_ERROR_CALLBACK(error_kind, ...)
#else
    #define GIF_ERROR_CALLBACK gif_error_callback_impl
#endif

void file_write_all(FILE *file, u8 const *source, isize size) {
    isize bytes_written = (isize) fwrite(source, 1, (size_t) size, file);

    // > On success, fread() and fwrite() return the number of items read or written.
    // > This number equals the number of bytes transferred only when size is 1.
    // > If an error occurs, or the end of the file is reached, the return value is a short item
    // > count (or zero).
    //
    // Which means that (size != bytes_written) indicates an error while writing to the file.

    if (size != bytes_written) {
        GIF_ERROR_CALLBACK(GIF_OUTPUT_IMAGE_WRITE_ERROR);
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
        GIF_ERROR_CALLBACK(GIF_OUT_OF_MEMORY);
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

void *alloc_resize(void *memory, isize old_size, isize new_size) {
    if (memory == NULL) {
        assert(old_size == 0);
        return alloc(new_size);
    }

    if (new_size == 0) {
        dealloc(memory, old_size);
        return NULL;
    }

    void *resized_memory = realloc(memory, (size_t) new_size);
    if (resized_memory == NULL) {
        GIF_ERROR_CALLBACK(GIF_OUT_OF_MEMORY);
        abort();
    }

    return resized_memory;
}

#define GIF_MAX_WIDTH 0xffff
#define GIF_MAX_HEIGHT 0xffff

typedef u8 ColorIndex;

typedef struct {
    ColorIndex *indices;
    isize count;
} ColorSequence;

u32 color_sequence_hash(ColorSequence sequence) {
    u32 hash = 5381;

    ColorIndex *sequence_iter = sequence.indices;
    while (sequence_iter < sequence.indices + sequence.count) {
        hash = ((hash << 5) + hash) + *sequence_iter;

        sequence_iter += 1;
    }

    return hash;
}

bool color_sequence_equals(ColorSequence left, ColorSequence right) {
    assert(left.count != 0 || right.count != 0);

    if (left.count != right.count) {
        return false;
    }

    return memcmp(left.indices, right.indices, (size_t) (left.count * sizeof(ColorIndex))) == 0;
}

typedef struct {
    ColorIndex *indices;
    isize count;
    isize capacity;
} ColorArray;

void color_array_push(ColorArray *colors, ColorIndex index) {
    if (colors->count == colors->capacity) {
        isize new_capacity = colors->capacity * 2;
        if (new_capacity < 16) {
            new_capacity = 16;
        }

        colors->indices = alloc_resize(
            colors->indices,
            colors->capacity * sizeof(ColorIndex),
            new_capacity * sizeof(ColorIndex)
        );
        colors->capacity = new_capacity;
    }

    colors->indices[colors->count] = index;
    colors->count += 1;
}

void color_array_destroy(ColorArray *colors) {
    dealloc(colors->indices, colors->capacity * sizeof(ColorIndex));
    colors->indices = NULL;
    colors->count = 0;
    colors->capacity = 0;
}

typedef isize LzwCode;
#define MAX_LZW_CODE 4095
#define LZW_CODE_MAX_BIT_LENGTH 12

// Key = ColorSequence
// Value = LzwCode
// An empty color sequence indicates an empty bucket.
typedef struct {
    ColorSequence sequence;
    LzwCode code;
} LzwTableBucket;

typedef struct {
    isize sequences_count;
    LzwTableBucket buckets[2 * (MAX_LZW_CODE + 1)];
} LzwTable;

LzwTable *lzw_table_create(void) {
    LzwTable *lzw_table = alloc(sizeof(LzwTable));
    memset(lzw_table, 0, (size_t) sizeof(LzwTable));

    return lzw_table;
}

void lzw_table_insert(LzwTable *lzw_table, ColorSequence new_sequence, LzwCode new_code) {
    assert(new_sequence.count > 0);
    assert(lzw_table->sequences_count < countof(lzw_table->buckets));

    isize bucket_index = color_sequence_hash(new_sequence) % countof(lzw_table->buckets);
    while (lzw_table->buckets[bucket_index].sequence.count > 0) {
        assert(!color_sequence_equals(new_sequence, lzw_table->buckets[bucket_index].sequence));

        bucket_index = (bucket_index + 1) % countof(lzw_table->buckets);
    }

    lzw_table->buckets[bucket_index].sequence = new_sequence;
    lzw_table->buckets[bucket_index].code = new_code;

    lzw_table->sequences_count += 1;
}

bool lzw_table_get_code(LzwTable *lzw_table, ColorSequence needle_sequence, LzwCode *code) {
    isize bucket_index = color_sequence_hash(needle_sequence) % countof(lzw_table->buckets);
    while (
        lzw_table->buckets[bucket_index].sequence.count > 0 &&
        !color_sequence_equals(needle_sequence, lzw_table->buckets[bucket_index].sequence)
    ) {
        bucket_index = (bucket_index + 1) % countof(lzw_table->buckets);
    }

    if (lzw_table->buckets[bucket_index].sequence.count > 0) {
        *code = lzw_table->buckets[bucket_index].code;
        return true;
    } else {
        return false;
    }
}

void lzw_table_destroy(LzwTable *lzw_table) {
    dealloc(lzw_table, sizeof(LzwTable));
}

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

#define RED_SHIFT 16
#define GREEN_SHIFT 8
#define BLUE_SHIFT 0

typedef u32 RGB;

RGB rgb(u8 red, u8 green, u8 blue) {
    return (RGB) red << RED_SHIFT | (RGB) green << GREEN_SHIFT | (RGB) blue << BLUE_SHIFT;
}

typedef struct {
    f32 x;
    f32 y;
    f32 z;
} f32x3;

f32x3 rgb_to_f32x3(RGB color) {
    f32x3 vector;
    vector.x = (f32) ((color >> RED_SHIFT) & 0xff) / 255.0F;
    vector.y = (f32) ((color >> GREEN_SHIFT) & 0xff) / 255.0F;
    vector.z = (f32) ((color >> BLUE_SHIFT) & 0xff) / 255.0F;

    return vector;
}

RGB f32x3_to_rgb(f32x3 vector) {
    RGB color = 0x000000;
    color |= (RGB) ((u8) (vector.x * 255.0F)) << RED_SHIFT;
    color |= (RGB) ((u8) (vector.y * 255.0F)) << GREEN_SHIFT;
    color |= (RGB) ((u8) (vector.z * 255.0F)) << BLUE_SHIFT;

    return color;
}

f32x3 f32x3_pow(f32x3 vector, f32 exponent) {
    f32x3 result;
    result.x = powf(vector.x, exponent);
    result.y = powf(vector.y, exponent);
    result.z = powf(vector.z, exponent);

    return result;
}

f32x3 f32x3_sub(f32x3 left, f32x3 right) {
    f32x3 result;
    result.x = left.x - right.x;
    result.y = left.y - right.y;
    result.z = left.z - right.z;

    return result;
}

f32 f32x3_dot(f32x3 left, f32x3 right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

#define GAMMA 2.2F

#define GIF_MAX_COLORS 256

typedef struct ColorTreeNode ColorTreeNode;

typedef struct {
    f32x3 color;
    ColorIndex color_index;
} IndexedColor;

struct ColorTreeNode {
    f32 median;
    ColorTreeNode *left;
    ColorTreeNode *right;
    IndexedColor value;
};

int indexed_color_compare_by_x(void const *left_void, void const *right_void) {
    f32x3 left = ((IndexedColor *) left_void)->color;
    f32x3 right = ((IndexedColor *) right_void)->color;

    return (int) copysignf(1.0F, left.x - right.x);
}

int indexed_color_compare_by_y(void const *left_void, void const *right_void) {
    f32x3 left = ((IndexedColor *) left_void)->color;
    f32x3 right = ((IndexedColor *) right_void)->color;

    return (int) copysignf(1.0F, left.y - right.y);
}

int indexed_color_compare_by_z(void const *left_void, void const *right_void) {
    f32x3 left = ((IndexedColor *) left_void)->color;
    f32x3 right = ((IndexedColor *) right_void)->color;

    return (int) copysignf(1.0F, left.z - right.z);
}

ColorTreeNode *color_tree_construct(IndexedColor *colors, isize colors_count, isize depth) {
    assert(colors_count != 0);

    ColorTreeNode *node = alloc(sizeof(ColorTreeNode));

    if (colors_count > 1) {
        isize dimension = depth % 3;

        int (*comparator)(void const *, void const *);
        switch (dimension) {
            case 0: comparator = indexed_color_compare_by_x; break;
            case 1: comparator = indexed_color_compare_by_y; break;
            case 2: comparator = indexed_color_compare_by_z; break;
            default: assert(false);
        }
        qsort(colors, (size_t) colors_count, sizeof(IndexedColor), comparator);

        IndexedColor median_color = colors[colors_count / 2];
        f32 median_color_array[3] = {
            median_color.color.x,
            median_color.color.y,
            median_color.color.z,
        };

        // The median value goes into the right subtree.
        node->median = median_color_array[dimension];

        node->left = color_tree_construct(
            colors,
            colors_count / 2,
            depth + 1
        );
        node->right = color_tree_construct(
            colors + colors_count / 2,
            colors_count - colors_count / 2,
            depth + 1
        );
    } else {
        // Both are NULL -> leaf node.
        node->left = NULL;
        node->right = NULL;

        node->value = colors[0];
    }

    return node;
}

void color_tree_get_nearest(
    ColorTreeNode *node,
    f32x3 search_color,
    isize depth,
    f32 *min_distance,
    ColorIndex *nearest_color_index
) {
    f32 search_color_array[3] = {search_color.x, search_color.y, search_color.z};

    if (node->left != NULL) {
        isize dimension = depth % 3;
        f32 dimension_value = search_color_array[dimension];

        ColorTreeNode *subtree;
        ColorTreeNode *other_subtree;
        if (dimension_value < node->median) {
            subtree = node->left;
            other_subtree = node->right;
        } else {
            subtree = node->right;
            other_subtree = node->left;
        }

        color_tree_get_nearest(
            subtree,
            search_color,
            depth + 1,
            min_distance,
            nearest_color_index
        );

        if (fabsf(node->median - dimension_value) < *min_distance) {
            color_tree_get_nearest(
                other_subtree,
                search_color,
                depth + 1,
                min_distance,
                nearest_color_index
            );
        }
    } else {
        // Leaf node
        // f32x3 distance_vector = f32x3_sub(node->value.color, search_color);
        // f32 current_distance = sqrtf(f32x3_dot(distance_vector, distance_vector));

        f32 current_distance = sqrtf(
            (node->value.color.x - search_color.x) * (node->value.color.x - search_color.x) +
            (node->value.color.y - search_color.y) * (node->value.color.y - search_color.y) +
            (node->value.color.z - search_color.z) * (node->value.color.z - search_color.z)
        );

        if (current_distance < *min_distance) {
            *min_distance = current_distance;
            *nearest_color_index = node->value.color_index;
        }
    }
}

typedef struct {
    isize count;
    RGB srgb_colors[GIF_MAX_COLORS];
    f32x3 rgb_colors[GIF_MAX_COLORS];
    ColorTreeNode *root;
} ColorTable;

void color_table_from_array(RGB *colors, isize colors_count, ColorTable *color_table) {
    assert(colors_count <= GIF_MAX_COLORS);

    color_table->count = colors_count;

    for (isize color_index = 0; color_index < colors_count; color_index += 1) {
        f32x3 rgb = rgb_to_f32x3(colors[color_index]);
        color_table->rgb_colors[color_index] = f32x3_pow(rgb, GAMMA);
        color_table->srgb_colors[color_index] = colors[color_index];
    }

    IndexedColor indexed_colors[GIF_MAX_COLORS];
    for (isize color_index = 0; color_index < colors_count; color_index += 1) {
        indexed_colors[color_index].color_index = (ColorIndex) color_index;
        indexed_colors[color_index].color = color_table->rgb_colors[color_index];
    }

    color_table->root = color_tree_construct(indexed_colors, color_table->count, 0);
}

void color_table_black_and_white(ColorTable *color_table) {
    RGB colors[2];
    colors[0] = rgb(0x00, 0x00, 0x00);
    colors[1] = rgb(0xff, 0xff, 0xff);

    color_table_from_array(colors, countof(colors), color_table);
}

void color_table_monochrome(ColorTable *color_table) {
    RGB colors[256];
    for (int component = 0x00; component <= 0xff; component += 1) {
        colors[component] = rgb((u8) component, (u8) component, (u8) component);
    }

    color_table_from_array(colors, countof(colors), color_table);
}

void color_table_web_safe(ColorTable *color_table) {
    RGB colors[216];
    isize count = 0;
    for (int red = 0x00; red <= 0xff; red += 0x33) {
        for (int green = 0x00; green <= 0xff; green += 0x33) {
            for (int blue = 0x00; blue <= 0xff; blue += 0x33) {
                assert(count < countof(colors));
                colors[count] = rgb((u8) red, (u8) green, (u8) blue);

                count += 1;
            }
        }
    }

    color_table_from_array(colors, countof(colors), color_table);
}

ColorIndex color_table_get_index(ColorTable *color_table, RGB needle) {
    f32x3 needle_vector = rgb_to_f32x3(needle);
    needle_vector = f32x3_pow(needle_vector, GAMMA);

    ColorIndex nearest_index;
    f32 min_distance = INFINITY;
    color_tree_get_nearest(color_table->root, needle_vector, 0, &min_distance, &nearest_index);

    return nearest_index;
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

    hash = ((hash << 5) + hash) + ((value >> RED_SHIFT) & 0xff);
    hash = ((hash << 5) + hash) + ((value >> GREEN_SHIFT) & 0xff);
    hash = ((hash << 5) + hash) + ((value >> BLUE_SHIFT) & 0xff);

    return hash;
}

void color_set_insert(ColorSet *colors, RGB new_color) {
    f32 fill_factor;
    if (colors->buckets_count == 0) {
        fill_factor = 1.0F;
    } else {
        fill_factor = (f32) colors->colors_count / (f32) colors->buckets_count;
    }

    if (fill_factor > 0.75F) {
        isize new_buckets_count = colors->buckets_count * 2;
        if (new_buckets_count < 16) {
            new_buckets_count = 16;
        }

        ColorSetBucket *new_buckets = alloc(new_buckets_count * sizeof(ColorSetBucket));
        memset(new_buckets, 0, (size_t) (new_buckets_count * sizeof(ColorSetBucket)));

        ColorSet rehashed_colors = {
            .buckets = new_buckets,
            .buckets_count = new_buckets_count,
            .colors_count = 0,
        };

        ColorSetBucket *old_buckets_iter = colors->buckets;
        while (old_buckets_iter < colors->buckets + colors->buckets_count) {
            if (old_buckets_iter->has_value) {
                color_set_insert(&rehashed_colors, old_buckets_iter->color);
            }
            old_buckets_iter += 1;
        }

        dealloc(colors->buckets, colors->buckets_count * sizeof(ColorSetBucket));
        colors->buckets = rehashed_colors.buckets;
        colors->buckets_count = rehashed_colors.buckets_count;
        colors->colors_count = rehashed_colors.colors_count;
    }

    isize bucket_index = rgb_hash(new_color) % colors->buckets_count;
    while (
        colors->buckets[bucket_index].has_value &&
        colors->buckets[bucket_index].color != new_color
    ) {
        bucket_index = (bucket_index + 1) % colors->buckets_count;
    }

    if (!colors->buckets[bucket_index].has_value) {
        colors->buckets[bucket_index].color = new_color;
        colors->buckets[bucket_index].has_value = true;
        colors->colors_count += 1;
    }
}

void color_set_destroy(ColorSet *colors) {
    dealloc(colors->buckets, colors->buckets_count * sizeof(ColorSetBucket));
    colors->buckets = NULL;
    colors->buckets_count = 0;
    colors->colors_count = 0;
}

int rgb_compare_red_only(void const *left_ptr, void const *right_ptr) {
    RGB left = *((RGB *) left_ptr);
    RGB right = *((RGB *) right_ptr);

    return (int) ((left >> RED_SHIFT) & 0xff) - (int) ((right >> RED_SHIFT) & 0xff);
}

int rgb_compare_green_only(void const *left_ptr, void const *right_ptr) {
    RGB left = *((RGB *) left_ptr);
    RGB right = *((RGB *) right_ptr);

    return (int) ((left >> GREEN_SHIFT) & 0xff) - (int) ((right >> GREEN_SHIFT) & 0xff);
}

int rgb_compare_blue_only(void const *left_ptr, void const *right_ptr) {
    RGB left = *((RGB *) left_ptr);
    RGB right = *((RGB *) right_ptr);

    return (int) ((left >> BLUE_SHIFT) & 0xff) - (int) ((right >> BLUE_SHIFT) & 0xff);
}

void color_table_median_cut(Image *image, isize target_colors_count, ColorTable *color_table) {
    assert(target_colors_count <= GIF_MAX_COLORS);
    RGB quantized_colors[GIF_MAX_COLORS];

    RGB *colors;
    RGB *colors_end;
    {
        ColorSet color_set = {NULL};
        u8 *image_data_iter = image->data;
        while (image_data_iter < image->data + image_size(image)) {
            RGB color = rgb(
                image_data_iter[RED_INDEX],
                image_data_iter[GREEN_INDEX],
                image_data_iter[BLUE_INDEX]
            );
            color_set_insert(&color_set, color);

            image_data_iter += image->bytes_per_pixel;
        }

        colors = alloc(color_set.colors_count * sizeof(RGB));
        colors_end = colors + color_set.colors_count;

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

        color_set_destroy(&color_set);
    }

    typedef struct {
        RGB *begin;
        RGB *end;
    } Segment;

    typedef struct {
        isize first_index;
        isize last_index;
        isize count;
        Segment segments[GIF_MAX_COLORS];
    } SegmentsQueue;

    SegmentsQueue queue;
    queue.first_index = 0;
    queue.last_index = 0;
    queue.count = 1;
    queue.segments[queue.first_index] = (Segment) {colors, colors_end};

    while (queue.count < target_colors_count) {
        assert(queue.count > 0);
        Segment current_segment = queue.segments[queue.first_index];

        // The image probably has less colors then we are targeting. Break out of the loop because
        // segments are naturally sorted by their size in decreasing order, so the rest of the
        // segments are going to be of size 1 as well.
        if (current_segment.end - current_segment.begin == 1) {
            break;
        }

        queue.first_index = (queue.first_index + 1) % countof(queue.segments);
        queue.count -= 1;

        u8 min_red = 0xff, min_green = 0xff, min_blue = 0xff;
        u8 max_red = 0x00, max_green = 0x00, max_blue = 0x00;
        {
            RGB *colors_iter = current_segment.begin;
            while (colors_iter != current_segment.end) {
                u8 red = (*colors_iter >> RED_SHIFT) & 0xff;
                min_red = u8_min(min_red, red);
                max_red = u8_max(max_red, red);

                u8 green = (*colors_iter >> GREEN_SHIFT) & 0xff;
                min_green = u8_min(min_green, green);
                max_green = u8_max(max_green, green);

                u8 blue = (*colors_iter >> BLUE_SHIFT) & 0xff;
                min_blue = u8_min(min_blue, blue);
                max_blue = u8_max(max_blue, blue);

                colors_iter += 1;
            }
        }

        int(*comparator)(void const *, void const *);
        if (
            max_red - min_red >= max_green - min_green &&
            max_red - min_red >= max_blue - min_blue
        ) {
            comparator = rgb_compare_red_only;
        } else if (
            max_green - min_green >= max_red - min_red &&
            max_green - min_green >= max_blue - min_blue
        ) {
            comparator = rgb_compare_green_only;
        } else {
            comparator = rgb_compare_blue_only;
        }
        qsort(
            current_segment.begin,
            (size_t) (current_segment.end - current_segment.begin),
            sizeof(RGB),
            comparator
        );

        queue.segments[(queue.last_index + 1) % countof(queue.segments)] = (Segment) {
            .begin = current_segment.begin,
            .end = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
        };
        queue.segments[(queue.last_index + 2) % countof(queue.segments)] = (Segment) {
            .begin = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
            .end = current_segment.end,
        };
        queue.last_index = (queue.last_index + 2) % countof(queue.segments);
        queue.count += 2;
    }

    RGB *quantized_colors_iter = quantized_colors;
    isize segment_index = queue.first_index;

    while (true) {
        u32 red_sum = 0;
        u32 green_sum = 0;
        u32 blue_sum = 0;

        RGB *segment_iter = queue.segments[segment_index].begin;
        while (segment_iter < queue.segments[segment_index].end) {
            red_sum += (*segment_iter >> RED_SHIFT) & 0xff;
            green_sum += (*segment_iter >> GREEN_SHIFT) & 0xff;
            blue_sum += (*segment_iter >> BLUE_SHIFT) & 0xff;

            segment_iter += 1;
        }

        // FIXME: Take gamma correction into account, when averaging the colors?

        // Should fit into a u32, because there are only 2^24 colors total.
        u32 colors_count =
            (u32) (queue.segments[segment_index].end - queue.segments[segment_index].begin);
        RGB average = rgb(
            (u8) (red_sum / colors_count),
            (u8) (green_sum / colors_count),
            (u8) (blue_sum / colors_count)
        );
        *quantized_colors_iter = average;

        if (segment_index == queue.last_index) {
            break;
        }
        quantized_colors_iter += 1;
        segment_index = (segment_index + 1) % countof(queue.segments);
    }

    dealloc(colors, (colors_end - colors) * sizeof(RGB));

    color_table_from_array(quantized_colors, queue.count, color_table);
}

int main(void) {
    Image image;
    image_load("input.jpg", &image);

    // TODO: Instead of quantizing the image on the go while encoding the image, convert the whole
    // image to color indices before starting the encoding process.

    ColorTable color_table;
    color_table_median_cut(&image, GIF_MAX_COLORS, &color_table);

    if (image.width > GIF_MAX_WIDTH || image.height > GIF_MAX_HEIGHT) {
        GIF_ERROR_CALLBACK(GIF_ERROR_INPUT_IMAGE_TOO_BIG, image.width, image.height);
        abort();
    }

    FILE *output_file = fopen("out.gif", "wb");
    if (output_file == NULL) {
        GIF_ERROR_CALLBACK(GIF_OUTPUT_IMAGE_WRITE_ERROR);
        abort();
    }


    // === GIF header ===

    // 6 bytes          for the "GIF89a" header
    // 7 bytes          for a logical screen descriptor
    // 256 * 3 bytes    for a max size global color table
    u8 header[6 + 7 + GIF_MAX_COLORS * 3];

    u8 *header_iter = header;

    // GIF89a
    header_iter = u8_write('G', header_iter);
    header_iter = u8_write('I', header_iter);
    header_iter = u8_write('F', header_iter);
    header_iter = u8_write('8', header_iter);
    header_iter = u8_write('9', header_iter);
    header_iter = u8_write('a', header_iter);

    // Logical screen descriptor
    header_iter = u16_write_le((u16) image.width, header_iter);
    header_iter = u16_write_le((u16) image.height, header_iter);

    if (color_table.count > 256) {
        GIF_ERROR_CALLBACK(GIF_ERROR_COLOR_TABLE_OVERFLOW);
        abort();
    }

    i32 color_table_size_log2_ceil = i32_log2_ceil((i32) color_table.count);
    if (color_table_size_log2_ceil < 1) {
        color_table_size_log2_ceil = 1;
    }
    assert(1 <= color_table_size_log2_ceil && color_table_size_log2_ceil <= 8);

    header_iter = u8_write(0x80 | (u8) (color_table_size_log2_ceil - 1), header_iter);

    header_iter = u8_write(0, header_iter);
    header_iter = u8_write(0, header_iter);

    // Global color table
    {
        RGB *color_iter = color_table.srgb_colors;

        for (isize i = 0; i < (1 << color_table_size_log2_ceil); i += 1) {
            if (color_iter < color_table.srgb_colors + color_table.count) {
                header_iter = u8_write((*color_iter >> RED_SHIFT) & 0xff, header_iter);
                header_iter = u8_write((*color_iter >> GREEN_SHIFT) & 0xff, header_iter);
                header_iter = u8_write((*color_iter >> BLUE_SHIFT) & 0xff, header_iter);

                color_iter += 1;
            } else {
                header_iter = u8_write(0x00, header_iter);
                header_iter = u8_write(0x00, header_iter);
                header_iter = u8_write(0x00, header_iter);
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
    image_header_iter = u8_write(0x21, image_header_iter);
    image_header_iter = u8_write(0xf9, image_header_iter);
    image_header_iter = u8_write(4, image_header_iter);
    image_header_iter = u8_write(0x08, image_header_iter);
    image_header_iter = u16_write_le(100, image_header_iter);
    image_header_iter = u8_write(0, image_header_iter);
    image_header_iter = u8_write(0, image_header_iter);

    // Image descriptor
    image_header_iter = u8_write(0x2c, image_header_iter);
    image_header_iter = u16_write_le(0, image_header_iter);
    image_header_iter = u16_write_le(0, image_header_iter);
    image_header_iter = u16_write_le((u16) image.width, image_header_iter);
    image_header_iter = u16_write_le((u16) image.height, image_header_iter);
    image_header_iter = u8_write(0x00, image_header_iter);

    // Minimum LZW code bit length
    isize lzw_code_min_bit_length = color_table_size_log2_ceil;
    *(image_header_iter++) = (u8) lzw_code_min_bit_length;

    file_write_all(output_file, image_header, image_header_iter - image_header);


    // === Image data ===

    LzwCodesWriter lzw_writer;
    lzw_writer_create(output_file, &lzw_writer);

    LzwCode lzw_clear_code = 1 << lzw_code_min_bit_length;
    LzwCode lzw_end_code = lzw_clear_code + 1;

    isize lzw_code_bit_length = lzw_code_min_bit_length + 1;
    isize next_lzw_code = lzw_clear_code + 2;

    // Output a clear code because the spec says so:
    // > Encoders should output a Clear code as the first code of each image data stream.
    lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);

    LzwTable *lzw_table = lzw_table_create();

    for (isize color_index = 0; color_index < (1 << color_table_size_log2_ceil); color_index += 1) {
        ColorSequence single_color_sequence = {
            .count = 1,
            .indices = alloc(sizeof(ColorIndex))
        };
        single_color_sequence.indices[0] = (ColorIndex) color_index;

        lzw_table_insert(lzw_table, single_color_sequence, color_index);
    }

    u8 *image_iter = image.data;
    u8 *image_end = image.data + image_size(&image);

    ColorArray current_sequence = {NULL};

    assert(image_end - image_iter > 0);

    {
        u8 red = image_iter[RED_INDEX];
        u8 green = image_iter[GREEN_INDEX];
        u8 blue = image_iter[BLUE_INDEX];
        image_iter += image.bytes_per_pixel;
        color_array_push(&current_sequence, color_table_get_index(&color_table, rgb(red, green, blue)));
    }

    while (image_iter < image_end || current_sequence.count > 0) {
        bool nonexistent_sequence_found = false;
        LzwCode longest_sequence_lzw_code = 0;

        u8 *image_iter_rewind = image_iter;

        while (true) {
            // current_sequence is never empty at this point
            assert(current_sequence.count > 0);

            bool sequence_already_exists = lzw_table_get_code(
                lzw_table,
                (ColorSequence) {current_sequence.indices, current_sequence.count},
                &longest_sequence_lzw_code
            );

            if (!sequence_already_exists) {
                nonexistent_sequence_found = true;
                break;
            }

            if (image_iter < image_end) {
                u8 red = image_iter[RED_INDEX];
                u8 green = image_iter[GREEN_INDEX];
                u8 blue = image_iter[BLUE_INDEX];
                image_iter += image.bytes_per_pixel;

                color_array_push(
                    &current_sequence,
                    color_table_get_index(&color_table, rgb(red, green, blue))
                );
            } else {
                break;
            }
        }

        // > When the table is full, the encoder can chose to use the table as is, making no changes
        // > to it until the encoder chooses to clear it. The encoder during this time sends out
        // > codes that are of the maximum Code Size.
        //
        // Which means that you don't have to necessarily emit the clear code here.

        if (nonexistent_sequence_found && next_lzw_code > MAX_LZW_CODE) {
            lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);

            image_iter = image_iter_rewind;
            current_sequence.count = 1;

            LzwTableBucket *lzw_table_iter = lzw_table->buckets;
            while (lzw_table_iter < lzw_table->buckets + countof(lzw_table->buckets)) {
                if (lzw_table_iter->sequence.count > 1) {
                    lzw_table->sequences_count -= 1;
                    dealloc(
                        lzw_table_iter->sequence.indices,
                        lzw_table_iter->sequence.count * sizeof(ColorIndex)
                    );
                    lzw_table_iter->sequence.count = 0;
                    lzw_table_iter->sequence.indices = NULL;
                }

                lzw_table_iter += 1;
            }

            lzw_code_bit_length = lzw_code_min_bit_length + 1;
            next_lzw_code = lzw_clear_code + 2;

            continue;
        }

        lzw_write_code(&lzw_writer, longest_sequence_lzw_code, lzw_code_bit_length);

        if (nonexistent_sequence_found) {
            // In here current_sequence length is at least 2 because all possible sequences of
            // length 1 are guaranteed to exist in the table.


            if (next_lzw_code <= MAX_LZW_CODE) {
                ColorSequence new_sequence = {
                    .count = current_sequence.count,
                    .indices = alloc(current_sequence.count * sizeof(ColorIndex)),
                };
                memmove(
                    new_sequence.indices,
                    current_sequence.indices,
                    (size_t) (current_sequence.count * sizeof(ColorIndex))
                );

                // Increase the LZW code length here because the spec says so:
                // > Whenever the LZW code value would exceed the current code length,
                // > the code length is increased by one. The packing/unpacking of these
                // > codes must then be altered to reflect the new code length.

                if ((next_lzw_code & (next_lzw_code - 1)) == 0) {
                    lzw_code_bit_length += 1;
                }

                lzw_table_insert(lzw_table, new_sequence, next_lzw_code);
                next_lzw_code += 1;
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
