#include <stdlib.h>     // malloc
#include <stddef.h>     // size_t, NULL
#include <stdio.h>      // fopen, fwrite, fclose
#include <string.h>     // memset, memcpy
#include <math.h>       // sqrtf, roundf
#include <stdbool.h>    // bool, true, false
#include <time.h>       // time
#include <stdint.h>     // uint32_t, uint64_t
#include <assert.h>     // assert

#include "common.h"
#include "gif_encoder.h"

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

typedef struct {
    isize x;
    isize y;
} isizex2;

typedef struct {
    GifColorIndex *pixels;
    isize width;
    isize height;
} Image;

void hexagon_draw(
    Image const *image,
    isizex2 center,
    isize side,
    GifColorIndex color,
    bool fill
) {
    isize inner_radius = (isize)roundf((f32)side * sqrtf(3.0F) / 2.0F);

    {
        isize x = center.x - inner_radius;
        isize y = center.y + side / 2;

        for (isize i = 0; i < side; i += 1) {
            if (fill) {
                for (isize j = x; j <= x + 2 * inner_radius; j += 1) {
                    image->pixels[y * image->width + j] = color;
                }
            } else {
                image->pixels[y * image->width + x] = color;
                image->pixels[y * image->width + x + 2 * inner_radius] = color;
            }
            y -= 1;
        }
    }

    isizex2 line_start = {center.x - inner_radius, center.y - side / 2};
    isizex2 line_end = {center.x, center.y - side};

    isize line_dx = line_end.x - line_start.x;
    isize line_dy = line_start.y - line_end.y;
    isize distance = 2 * line_dy - line_dx;
    isize x = line_start.x;
    isize y = line_start.y;
    while (true) {
        if (fill) {
            for (isize i = x; i <= line_start.x + 2 * inner_radius - (x - line_start.x); i += 1) {
                image->pixels[y * image->width + i] = color;
                image->pixels[(line_start.y + side + (line_start.y - y)) * image->width + i] = color;
            }
        } else {
            image->pixels[y * image->width + x] = color;
            image->pixels[(line_start.y + side + (line_start.y - y)) * image->width + x] = color;
            image->pixels[y * image->width + line_start.x + 2 * inner_radius - (x - line_start.x)] = color;
            image->pixels[(line_start.y + side + (line_start.y - y)) * image->width + line_start.x + 2 * inner_radius - (x - line_start.x)] = color;
        }

        if (x == line_end.x) {
            break;
        }

        if (distance >= 0) {
            y -= 1;
            distance += 2 * (line_dy - line_dx);
        } else {
            distance += 2 * line_dy;
        }

        x += 1;
    }
}

int main(void) {
    pcg32_random_t rng;
    rng.state = 0U;
    rng.inc = ((uint64_t)time(NULL) << 1u) | 1u;
    pcg32_random_r(&rng);
    rng.state += 0;
    pcg32_random_r(&rng);

    isize arena_capacity = 128 * 1024 * 1024;
    u8 *arena_memory = malloc((size_t)arena_capacity);
    GifArena arena = {arena_memory, arena_memory + arena_capacity};

    isize grid_width = 120;
    isize grid_height = 120;
    isize side = 7;

    bool *cells_ping = gif_arena_alloc(&arena, grid_width * grid_height * sizeof(bool));
    bool *cells_pong = gif_arena_alloc(&arena, grid_width * grid_height * sizeof(bool));

    isize inner_radius = (isize)roundf((f32)side * sqrtf(3.0F) / 2.0F);

    Image image = {
        .width = grid_width * 2 * inner_radius + inner_radius + 1,
        .height = grid_height * (side / 2 + side + 1) + (side - side / 2) + 1,
    };
    image.pixels = gif_arena_alloc(&arena, image.width * image.height);

    isize color_count = srgb_palette_monochrome(NULL);
    u8 *colors = gif_arena_alloc(&arena, color_count * 3);
    srgb_palette_monochrome(colors);

    FILE *output_file = fopen("autohex.gif", "wb");

    GifEncoder *encoder = gif_encoder_create(&arena);

    GifOutputBuffer out_buffer;
    gif_out_buffer_create(8192, &out_buffer, &arena);

#if 0
    for (isize i = 0; i < grid_width * grid_height; i += 1) {
        u32 random = pcg32_random_r(&rng);
        if (random < ((u32)1 << 31)) {
            cells_ping[i] = false;
        } else {
            cells_ping[i] = true;
        }
    }
#endif
    memset(cells_ping, 0, (size_t)(grid_width * grid_height * sizeof(bool)));
    cells_ping[(grid_height - 1) * grid_width + (grid_width - 1)] = true;

    gif_encoder_start(encoder, image.width, image.height, colors, color_count, &out_buffer);
    for (int frame = 0; frame < 250; frame += 1) {
        memset(image.pixels, 0, (size_t)(image.width * image.height));
        isize y = side;
        for (isize grid_y = 0; grid_y < grid_height; grid_y += 1) {
            isize x;
            if (grid_y % 2 == 0) {
                x = inner_radius;
            } else {
                x = 2 * inner_radius;
            }

            for (isize grid_x = 0; grid_x < grid_width; grid_x += 1) {
                bool fill = cells_ping[grid_y * grid_width + grid_x];

                hexagon_draw(&image, (isizex2){x, y}, side, 255, fill);

                x += 2 * inner_radius;
            }

            y += 3 * side / 2 + 1;
        }

        GifColorIndex *image_iter = image.pixels;
        GifColorIndex *image_end = image.pixels + image.width * image.height;

        gif_encoder_start_frame(encoder, NULL, 0, &out_buffer);
        if (gif_out_buffer_capacity_left(&out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
            fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
            gif_out_buffer_reset(&out_buffer);
        }

        while (image_iter != image_end) {
            image_iter += gif_encoder_feed_frame(
                encoder,
                image_iter,
                image_end - image_iter,
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

        // rule 110 (bottom half)

        memset(cells_pong, 0, (size_t)(grid_width * grid_height * sizeof(bool)));

        for (isize grid_y = grid_height / 2; grid_y < grid_height - 1; grid_y += 1) {
            memcpy(&cells_ping[grid_y * grid_width], &cells_ping[(grid_y + 1) * grid_width], (size_t)(grid_width * sizeof(bool)));
            memcpy(&cells_pong[grid_y * grid_width], &cells_ping[(grid_y + 1) * grid_width], (size_t)(grid_width * sizeof(bool)));
        }

        for (isize grid_x = 0; grid_x < grid_width; grid_x += 1) {
            bool n0 = cells_ping[(grid_height - 2) * grid_width + (grid_x - 1 + grid_width) % grid_width];
            bool n1 = cells_ping[(grid_height - 2) * grid_width + grid_x];
            bool n2 = cells_ping[(grid_height - 2) * grid_width + (grid_x + 1) % grid_width];

            u32 pattern = (n0 ? 1 << 2 : 0) | (n1 ? 1 << 1 : 0) | (n2 ? 1 : 0);
            assert(0 <= pattern && pattern <= 7);

            bool table[] = {
                [7 /* 111 */] = 0,
                [6 /* 110 */] = 1,
                [5 /* 101 */] = 1,
                [4 /* 100 */] = 0,
                [3 /* 011 */] = 1,
                [2 /* 010 */] = 1,
                [1 /* 001 */] = 1,
                [0 /* 000 */] = 0,
            };

            cells_ping[(grid_height - 1) * grid_width + grid_x] = table[pattern];
            cells_pong[(grid_height - 1) * grid_width + grid_x] = table[pattern];
        }

        // game of life (top half)

        for (isize grid_index = 0; grid_index <= grid_width * (grid_height / 2); grid_index += 1) {
            isize current_x = grid_index % grid_width;
            isize current_y = grid_index / grid_width;

            isize alive_neighbors = 0;

            isize dxs[] = {
                -1, // west
                1, // east
                current_y % 2 == 0 ? -1 : 0, // north-west
                current_y % 2 == 0 ? 0 : 1, // north-east
                current_y % 2 == 0 ? -1 : 0, // south-west
                current_y % 2 == 0 ? 0 : 1, // south-east
            };
            isize dys[] = {
                0, // west
                0, // east
                -1, // north-west
                -1, // north-east
                1, // south-west
                1, // south-east
            };

            for (isize i = 0; i < 6; i += 1) {
                // isize neighbor_x = (current_x + dxs[i] + grid_width) % grid_width;
                // isize neighbor_y = (current_y + dys[i] + grid_height) % grid_height;
                isize neighbor_x = (current_x + dxs[i] + grid_width) % grid_width;

                isize neighbor_y = current_y + dys[i];
                if (neighbor_y < 0) {
                    continue;
                }

                if (cells_ping[neighbor_y * grid_width + neighbor_x]) {
                    alive_neighbors += 1;
                }
            }

            if (cells_ping[current_y * grid_width + current_x]) {
                if (alive_neighbors == 2) {
                    cells_pong[current_y * grid_width + current_x] = true;
                } else {
                    cells_pong[current_y * grid_width + current_x] = false;
                }
            } else {
                if (alive_neighbors == 2) {
                    cells_pong[current_y * grid_width + current_x] = true;
                } else {
                    cells_pong[current_y * grid_width + current_x] = false;
                }
            }
        }
        {
            bool *swap = cells_ping;
            cells_ping = cells_pong;
            cells_pong = swap;
        }
    }
    gif_encoder_finish(encoder, &out_buffer);

    // Write any leftover buffered data into the file.
    fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
    fclose(output_file);

    return 0;
}
