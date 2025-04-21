#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "gif_encoder.h"

int main(void) {
    isize arena_capacity = 128 * 1024 * 1024;
    u8 *arena_memory = malloc((size_t)arena_capacity);
    GifArena arena = {arena_memory, arena_memory + arena_capacity};

    FILE *input_file = fopen("./debug/missing_last_lzw_code.gif", "rb");
    if (input_file == NULL) {
        return 1;
    }

    fseek(input_file, 0, SEEK_END);
    long input_file_size = ftell(input_file);
    fseek(input_file, 0, SEEK_SET);
    printf("input file size: %ld\n", input_file_size);

    u8 *gif = gif_arena_alloc(&arena, input_file_size);
    u8 *gif_end = gif + input_file_size;
    unsigned long bytes_read = fread(gif, 1, (size_t)input_file_size, input_file);
    assert((long)bytes_read == input_file_size);
    assert(memcmp(gif, "GIF89a", 6) == 0);

    FILE *output_file = fopen("./debug/missing_last_lzw_code_trimmed.gif", "wb");
    if (output_file == NULL) {
        return 1;
    }

    u8 *gif_iter = gif;
    u8 *gif_chunk_begin = gif;
    gif_iter += 6; // GIF89a
    gif_iter += 7; // logical screen descriptor
    gif_iter += 256 * 3; // global color table

    fwrite(gif_chunk_begin, 1, (size_t)(gif_iter - gif_chunk_begin), output_file);
    gif_chunk_begin = gif_iter;

    isize frame_index = 1;
    while (*gif_iter != 0x3b) {
        gif_iter += 8; // graphic control extension
        gif_iter += 10; // image descriptor
        gif_iter += 1; // min LZW code size

        while (*gif_iter != 0) {
            isize block_size = *gif_iter;
            gif_iter += 1;
            gif_iter += block_size;
        }
        gif_iter += 1;

        if (frame_index == 58) {
            fwrite(gif_chunk_begin, 1, (size_t)(gif_iter - gif_chunk_begin), output_file);
        }
        gif_chunk_begin = gif_iter;

        frame_index += 1;
    }

    gif_iter += 1;
    fwrite(gif_chunk_begin, 1, (size_t)(gif_iter - gif_chunk_begin), output_file);
    assert(gif_iter == gif_end);

    fclose(output_file);
    return 0;
}
