// $ cc -DNDEBUG -O2 main.c
//
// GIF references:
// - https://www.w3.org/Graphics/GIF/spec-gif89a.txt
// - https://en.wikipedia.org/wiki/GIF

// TODO: Come up with a GIF encoder interface and extract it into a separate file.

#include <stdio.h>      // FILE, fopen, fclose, fwrite, fprintf, stderr
#include <stdint.h>     // uint8_t, uint16_t, ptrdiff_t
#include <stddef.h>     // size_t, NULL
#include <stdlib.h>     // malloc, free, exit
#include <string.h>     // memmove, memcmp, memset
#include <stdbool.h>    // true, false, bool
#include <assert.h>     // assert

#define UNUSED(x) ((void) (x))

typedef uint8_t u8;
typedef uint16_t u16;
typedef ptrdiff_t isize;

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
        exit(1);
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
        exit(1);
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

#define IMAGE_BYTES_PER_PIXEL 4
#define RED_INDEX 2
#define GREEN_INDEX 1
#define BLUE_INDEX 0

int main(void) {
    isize image_width = 800;
    isize image_height = 600;
    if (image_width > GIF_MAX_WIDTH || image_height > GIF_MAX_HEIGHT) {
        // Image dimensions are too big for a GIF.
        exit(1);
    }

    u8 *image = alloc(image_width * image_height * IMAGE_BYTES_PER_PIXEL);
    /* Generate a blank white image */ {
        u8 *image_iter = image;
        for (isize y = 0; y < image_height; y += 1) {
            for (isize x = 0; x < image_width; x += 1) {
                image_iter[RED_INDEX] = 0xff;
                image_iter[GREEN_INDEX] = 0xff;
                image_iter[BLUE_INDEX] = 0xff;

                image_iter += IMAGE_BYTES_PER_PIXEL;
            }
        }
    }

    FILE *output_file = fopen("out.gif", "wb");
    if (output_file == NULL) {
        // Failed to open an output file.
        exit(1);
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
    header_iter = U8_WRITE(0x80, header_iter);
    header_iter = U8_WRITE(0, header_iter);
    header_iter = U8_WRITE(0, header_iter);

    // Global color table (black and white)
    header_iter = U8_WRITE(0x00, header_iter);
    header_iter = U8_WRITE(0x00, header_iter);
    header_iter = U8_WRITE(0x00, header_iter);
    header_iter = U8_WRITE(0xff, header_iter);
    header_iter = U8_WRITE(0xff, header_iter);
    header_iter = U8_WRITE(0xff, header_iter);

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
    isize lzw_code_min_bit_length = 2;
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

    ColorIndex black_index = 0;
    ColorSequence *black_color_sequence = alloc(sizeof(ColorSequence) + sizeof(ColorIndex));
    black_color_sequence->count = 1;
    black_color_sequence->indices[0] = black_index;
    lzw_table->sequences[black_index] = black_color_sequence;

    ColorIndex white_index = 1;
    ColorSequence *white_color_sequence = alloc(sizeof(ColorSequence) + sizeof(ColorIndex));
    white_color_sequence->count = 1;
    white_color_sequence->indices[0] = white_index;
    lzw_table->sequences[white_index] = white_color_sequence;

    // Create a dummy sequence for the unused LZW codes, as well as for the "clear" and the "end of
    // information" codes, so that we don't have to handle these codes differently from others when
    // looking up sequences in the LZW table.
    ColorSequence *dummy_sequence = alloc(sizeof(ColorSequence));
    dummy_sequence->count = 0;

    // TODO: Start from the "size of the color table" instead of a hard-coded value
    // (2 for the black-and-white table).
    for (LzwCode dummy_code = 2; dummy_code < next_lzw_code; dummy_code += 1) {
        lzw_table->sequences[dummy_code] = dummy_sequence;
    }

    u8 *image_iter = image;
    u8 *image_end = image + image_width * image_height * IMAGE_BYTES_PER_PIXEL;

    ColorArray current_sequence = {NULL};

    assert(image_end - image_iter > 0);

    // TODO: Abstract away looking up color indices from the color table.

    u8 red = image_iter[RED_INDEX];
    u8 green = image_iter[GREEN_INDEX];
    u8 blue = image_iter[BLUE_INDEX];
    image_iter += IMAGE_BYTES_PER_PIXEL;

    if (red == 0x00 && green == 0x00 && blue == 0x00) {
        color_array_push(&current_sequence, black_index);
    } else if (red == 0xff && green == 0xff && blue == 0xff) {
        color_array_push(&current_sequence, white_index);
    } else {
        // TODO: Search for the closest available color? Preprocess the image so that it only
        // contains colors from the color table?
        exit(1);
    }

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
                // TODO: Code duplication.

                u8 blue = *(image_iter + 0);
                u8 green = *(image_iter + 1);
                u8 red = *(image_iter + 2);
                image_iter += IMAGE_BYTES_PER_PIXEL;

                if (red == 0x00 && green == 0x00 && blue == 0x00) {
                    color_array_push(&current_sequence, black_index);
                } else if (red == 0xff && green == 0xff && blue == 0xff) {
                    color_array_push(&current_sequence, white_index);
                } else {
                    exit(1);
                }
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
            image_iter = image_iter_rewind;
            current_sequence.count = 1;

            for (LzwCode i = lzw_clear_code + 2; i < lzw_table->count; i += 1) {
                ColorSequence *sequence = lzw_table->sequences[i];
                dealloc(sequence, sizeof(ColorSequence) + sequence->count * sizeof(ColorIndex));
            }
            lzw_table->count = lzw_clear_code + 2;
            lzw_code_bit_length = lzw_code_min_bit_length + 1;

            lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);
            continue;
        }

        lzw_write_code(&lzw_writer, longest_sequence_lzw_code, lzw_code_bit_length);

        if (nonexistent_sequence_found) {
            // In here current_sequence length is at least 2 because all possible sequences of
            // length 1 are guaranteed to exist in the table.

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
