#include <stdio.h> // FILE, fopen, fclose, fwrite
#include <stdint.h> // uint8_t, uint16_t
#include <stddef.h> // size_t, NULL
#include <stdlib.h> // malloc, free, exit
#include <string.h> // memmove, memcmp
#include <stdbool.h> // true, false, bool
#include <assert.h> // assert

#define BYTES_PER_PIXEL 4

size_t u16_write_le(uint16_t value, uint8_t *dest) {
    *(dest + 0) = (uint8_t) (value & 0xff);
    *(dest + 1) = (uint8_t) (value >> 8);

    return 2;
}

#define MAX_LZW_CODE 4095

typedef uint8_t ColorIndex;

typedef struct {
    size_t count;
    ColorIndex indices[];
} ColorsSequence;

typedef struct {
    ColorsSequence *sequences[MAX_LZW_CODE + 1];
    size_t count;
} LzwTable;

typedef struct {
    ColorIndex *indices;
    size_t count;
    size_t capacity;
} DynamicColorsSequence;

void colors_sequece_push(DynamicColorsSequence *sequence, ColorIndex index) {
    if (sequence->count == sequence->capacity) {
        size_t new_capacity = sequence->capacity * 2;
        if (new_capacity == 0) {
            new_capacity = 16;
        }

        ColorIndex *new_indices = malloc(new_capacity * sizeof(ColorIndex));
        memmove(new_indices, sequence->indices, sequence->count);
        free(sequence->indices);
        sequence->indices = new_indices;
        sequence->capacity = new_capacity;
    }

    sequence->indices[sequence->count] = index;
    sequence->count += 1;
}

typedef struct {
    FILE *output_file;
    uint8_t block[256];
    uint8_t *block_iter;
    size_t bit_pos;
} LzwWriter;

void lzw_writer_create(FILE *output_file, LzwWriter *writer) {
    writer->output_file = output_file;
    // skip over the size byte (we fill it later in the lzw_flush function)
    writer->block_iter = writer->block + 1;
    writer->bit_pos = 0;

    memset(writer->block, 0, 256);
}

void lzw_flush(LzwWriter *writer) {
    // the block is empty no need to output it to the file
    if (writer->block_iter == writer->block + 1 && writer->bit_pos == 0) {
        return;
    }

    uint8_t block_size = writer->block_iter - (writer->block + 1);
    if (writer->bit_pos != 0) {
        block_size += 1;
    }
    writer->block[0] = block_size;

    size_t items_written = fwrite(writer->block, sizeof(uint8_t), block_size + 1, writer->output_file);
    assert(items_written == block_size + 1);

    // skip over the size byte (we fill it later in the lzw_flush function)
    writer->block_iter = writer->block + 1;
    writer->bit_pos = 0;

    memset(writer->block, 0, 256);
}

#define BLOCK_SIZE 256

void lzw_write_code(LzwWriter *writer, size_t code, size_t code_bit_length) {
    // printf("%zu\n", code);

    // size_t full_bytes_left = (writer->block + 2) - (writer->block_iter + 1);
    // size_t writer_bits_left = 8 - writer->bit_pos + 8 * full_bytes_left;
    // if (writer_bits_left < code_bit_length) {
    //     lzw_flush(writer);
    // }

    if (writer->block + BLOCK_SIZE == writer->block_iter) {
        lzw_flush(writer);
    }

    size_t code_bits_left = code_bit_length;
    size_t code_left = code;
    
    *writer->block_iter |= (uint8_t) ((code_left << writer->bit_pos) & 0xff);
    if (code_bits_left + writer->bit_pos < 8) {
        writer->bit_pos += code_bits_left;
        code_bits_left = 0;
    } else {
        code_left >>= (8 - writer->bit_pos);
        code_bits_left -= (8 - writer->bit_pos);
        writer->bit_pos = 0;
        writer->block_iter += 1;
        /**writer->block_iter = 0x00;*/
    }

    // executed a max of 2 times because in the worst case
    // we're gonna have to write 12 bits across 3 bytes
    while (code_bits_left > 0) {
        if (writer->block + BLOCK_SIZE == writer->block_iter) {
            lzw_flush(writer);
        }

        // writer->bit_pos is guaranteed to be 0
        *writer->block_iter |= (uint8_t) (code_left & 0xff);
        
        if (code_bits_left < 8) {
            writer->bit_pos = code_bits_left;
            code_bits_left = 0;
        } else {
            code_left >>= 8;
            code_bits_left -= 8;
            writer->bit_pos = 0;
            writer->block_iter += 1;
            /**writer->block_iter = 0x00;*/
        }
    }
}

int main(void) {
    uint16_t width = 4096;
    uint16_t height = 4096;
    uint8_t *image = malloc((size_t) width * (size_t) height * BYTES_PER_PIXEL);
    {
        uint8_t *image_iter = image;
        for (size_t y = 0; y < height; y += 1) {
            for (size_t x = 0; x < width; x += 1) {
                *(image_iter++) = 0xff; // blue
                *(image_iter++) = 0xff; // green
                *(image_iter++) = 0xff; // red
                *(image_iter++) = 0x00;
            }
        }
    }

    FILE *output_file = fopen("out.gif", "w");

    // 6 bytes for header
    // 7 bytes for logical screen descriptor
    // 256 * 3 bytes for max size global color table
    uint8_t header[6 + 7 + 256 * 3];
    uint8_t *header_iter = header;

    // header
    *(header_iter++) = 'G';
    *(header_iter++) = 'I';
    *(header_iter++) = 'F';
    *(header_iter++) = '8';
    *(header_iter++) = '9';
    *(header_iter++) = 'a';

    // logical screen descriptor
    header_iter += u16_write_le(width, header_iter);
    header_iter += u16_write_le(height, header_iter);
    *(header_iter++) = 0x80;
    *(header_iter++) = 0x00;
    *(header_iter++) = 0x00;

    // global color table (black and white)
    *(header_iter++) = 0x00;
    *(header_iter++) = 0x00;
    *(header_iter++) = 0x00;
    *(header_iter++) = 0xff;
    *(header_iter++) = 0xff;
    *(header_iter++) = 0xff;

    {
        size_t items_written = fwrite(header, sizeof(uint8_t), header_iter - header, output_file);
        assert(items_written == header_iter - header);
    }

    // Image
    uint8_t image_descriptor[10];
    uint8_t *image_descriptor_iter = image_descriptor;

    *(image_descriptor_iter++) = 0x2c;
    image_descriptor_iter += u16_write_le(0, image_descriptor_iter);
    image_descriptor_iter += u16_write_le(0, image_descriptor_iter);
    image_descriptor_iter += u16_write_le(width, image_descriptor_iter);
    image_descriptor_iter += u16_write_le(height, image_descriptor_iter);
    *(image_descriptor_iter++) = 0x0;

    {
        size_t items_written = fwrite(image_descriptor, sizeof(uint8_t), sizeof(image_descriptor), output_file);
        assert(items_written == sizeof(image_descriptor));
    }

    // Image data
    uint8_t min_lzw_code[] = {2};
    fwrite(min_lzw_code, sizeof(uint8_t), sizeof(min_lzw_code), output_file);

    uint8_t image_sub_block[256];
    uint8_t *image_sub_block_iter = image_sub_block;

    LzwWriter lzw_writer;
    lzw_writer_create(output_file, &lzw_writer);
    
    size_t lzw_code_bit_length = min_lzw_code[0] + 1;

    // 2 ^ min_lzw_code
    size_t lzw_clear_code = 4;
    size_t lzw_end_code = 5;

    // Output a clear code because the spec says so:
    // > Encoders should output a Clear code as the first code of each image data stream.
    lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);

    LzwTable *lzw_table = malloc(sizeof(LzwTable));

    ColorIndex black_index = 0;
    ColorsSequence *black_color = malloc(sizeof(ColorsSequence) + sizeof(ColorIndex));
    black_color->count = 1;
    black_color->indices[0] = 0;
    lzw_table->sequences[black_index] = black_color;

    ColorIndex white_index = 1;
    ColorsSequence *white_color = malloc(sizeof(ColorsSequence) + sizeof(ColorIndex));
    white_color->count = 1;
    white_color->indices[0] = 1;
    lzw_table->sequences[white_index] = white_color;

    ColorsSequence *dummy_sequence = malloc(sizeof(ColorsSequence));
    dummy_sequence->count = 0;
    lzw_table->sequences[2] = dummy_sequence;
    lzw_table->sequences[3] = dummy_sequence;
    lzw_table->sequences[4] = dummy_sequence;
    lzw_table->sequences[5] = dummy_sequence;

    lzw_table->count = 6;

    uint8_t *image_iter = image;
    uint8_t *image_end = image + (size_t) width * (size_t) height * BYTES_PER_PIXEL;

    DynamicColorsSequence current_sequence = {0};

    // TODO: assert that an image has at least a single pixel

    uint8_t blue = *(image_iter + 0);
    uint8_t green = *(image_iter + 1);
    uint8_t red = *(image_iter + 2);
    image_iter += BYTES_PER_PIXEL;

    if (red == 0x00 && green == 0x00 && blue == 0x00) {
        colors_sequece_push(&current_sequence, black_index);
    } else if (red == 0xff && green == 0xff && blue == 0xff) {
        colors_sequece_push(&current_sequence, white_index);
    } else {
        exit(1);
    }

    while (image_iter < image_end || current_sequence.count > 0) {
        bool nonexistent_sequence_found = false;
        size_t longest_sequence_lzw_code = 0;

        uint8_t *image_iter_rewind = image_iter;

        while (true) {
            bool sequence_already_exists = false;

            for (size_t lzw_code = 0; lzw_code < lzw_table->count; lzw_code += 1) {
                if (current_sequence.count != lzw_table->sequences[lzw_code]->count) {
                    continue;
                }

                // current_sequence always has length > 1, so memcmp is not going to return 0 for dummy sequences
                if (memcmp(current_sequence.indices, lzw_table->sequences[lzw_code]->indices, current_sequence.count) == 0) {
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
                uint8_t blue = *(image_iter + 0);
                uint8_t green = *(image_iter + 1);
                uint8_t red = *(image_iter + 2);
                image_iter += BYTES_PER_PIXEL;

                if (red == 0x00 && green == 0x00 && blue == 0x00) {
                    colors_sequece_push(&current_sequence, black_index);
                } else if (red == 0xff && green == 0xff && blue == 0xff) {
                    colors_sequece_push(&current_sequence, white_index);
                } else {
                    exit(1);
                }
            } else {
                break;
            }
        }

        // input color indices
        // 1 | 1 1 | 1 | 1 1 | 1 | 1 1 | 1 | 1 1 | 1 1 1
        //
        // encoded output
        // 1 6 4 1 6 4 1 6 4 1 6 7 5
        //
        // lzw max code = 7
        //
        // lzw table
        // 0 -> 0
        // 1 -> 1
        // 2 -> dummy
        // 3 -> dummy
        // 4 -> clear code
        // 5 -> end
        // 6 -> 1 1
        // 7 -> 1 1 1

        if (nonexistent_sequence_found && lzw_table->count == MAX_LZW_CODE + 1) {
            image_iter = image_iter_rewind;
            current_sequence.count = 1;

            for (size_t i = 6; i < lzw_table->count; i += 1) {
                free(lzw_table->sequences[i]);
            }
            lzw_table->count = 6;
            lzw_code_bit_length = 3;

            lzw_write_code(&lzw_writer, lzw_clear_code, lzw_code_bit_length);
            continue;
        }

        // Emit the longest_sequence_lzw_code to the output here
        lzw_write_code(&lzw_writer, longest_sequence_lzw_code, lzw_code_bit_length);

        if (nonexistent_sequence_found) {
            // in here current_sequence length is at least 2 because all possible sequences of length
            // 1 are guaranteed to exist in the table

            ColorsSequence *new_sequence = malloc(sizeof(ColorsSequence) + current_sequence.count * sizeof(ColorIndex));
            new_sequence->count = current_sequence.count;
            memmove(new_sequence->indices, current_sequence.indices, current_sequence.count * sizeof(ColorIndex));

            // for (size_t i = 0; i < new_sequence->count; i += 1) {
            //     printf("%d", new_sequence->indices[i]);
            //     if (new_sequence->indices[i] != 1) {
            //         exit(42);
            //     }
            // }
            // fprintf(stdout, "%zu\n", new_sequence->count);
            // fflush(stdout);

            lzw_table->sequences[lzw_table->count] = new_sequence;

            // increase the LZW code length here because the spec says so:
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
            // TODO: assert that we hit the end of the input data here
            current_sequence.count = 0;
        }
    }

    // Output an "end of information" code here
    lzw_write_code(&lzw_writer, lzw_end_code, lzw_code_bit_length);
    lzw_flush(&lzw_writer);

    uint8_t block_terminator[] = {0x00};
    fwrite(block_terminator, sizeof(uint8_t), sizeof(block_terminator), output_file);

    // Trailer
    uint8_t trailer[] = {0x3b};
    fwrite(trailer, sizeof(uint8_t), sizeof(trailer), output_file);

    fclose(output_file);
    free(image);

    for (size_t i = 6; i < lzw_table->count; i += 1) {
        free(lzw_table->sequences[i]);
    }

    free(lzw_table);

    free(white_color);
    free(black_color);
    free(dummy_sequence);

    free(current_sequence.indices);

    return 0;
}
