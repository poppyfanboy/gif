// TODO: Crop the GIF, so that the grid fully covers the canvas.
// TODO: Make the visualization animated: hexagons should fade in/out and expand/shrink smoothly.
// TODO: Paint cells which were previously alive with a different color slowly fading into black?
//
// I came up with / stole these while looking at https://github.com/bobombolo/hexomata screenshots:
// TODO: Try running simulation from a single alive cell. Random / rule 110 feeding is too chaotic.
// TODO: Killing the whole previous generation of cells is also a valid rule!
// TODO: Make the rules directional? (e.g. a cell is born only if it has west and east neighbors.)
// TODO: Take more neighbors into account? 18 instead of just 6? (+12 from the second layer around.)
// TODO: Make the rules randomly generated? Combine with the directional rules idea?
// TODO: Color cells depending on the number of neighbors or the rule which gave them birth?

#include <stdlib.h>     // malloc
#include <stddef.h>     // size_t, NULL
#include <stdio.h>      // fopen, fwrite, fclose, printf
#include <string.h>     // memset, memcpy
#include <math.h>       // sqrtf, roundf
#include <stdbool.h>    // bool, true, false
#include <time.h>       // time
#include <assert.h>     // assert

#include "common.h"
#include "gif_encoder.h"

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct {
    u64 state;
    u64 inc;
} PCG32;

u32 pcg32_random(PCG32 *rng) {
    u64 oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ull + (rng->inc | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    u32 xorshifted = (u32)(((oldstate >> 18u) ^ oldstate) >> 27u);
    u32 rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void pcg32_init(PCG32 *rng, u64 init_state) {
    rng->state = 0u;
    rng->inc = 1u;
    pcg32_random(rng);
    rng->state += init_state;
    pcg32_random(rng);
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

void draw_hexagon(Image *image, isizex2 center, isize side, GifColorIndex color) {
    // Even side length => left/right sides take up an odd number of pixels => pixel perfect center.
    assert(side % 2 == 0);

    isize radius = (isize)roundf((f32)side * sqrtf(3.0F) / 2.0F);

    // Draw west and east sides.
    for (isize y = center.y - side / 2; y <= center.y + side / 2; y += 1) {
        image->pixels[y * image->width + center.x - radius] = color;
        image->pixels[y * image->width + center.x + radius] = color;
    }

    // Draw the south-west side using Bresenham and mirror it both horizontally and vertically.
    isize dx = radius;
    isize dy = side / 2;
    isize error = 2 * dy - dx;

    isize y = 0;
    for (isize x = 0; x <= radius; x += 1) {
        image->pixels[(center.y - side / 2 - y) * image->width + (center.x - radius + x)] = color;
        image->pixels[(center.y - side / 2 - y) * image->width + (center.x + radius - x)] = color;
        image->pixels[(center.y + side / 2 + y) * image->width + (center.x - radius + x)] = color;
        image->pixels[(center.y + side / 2 + y) * image->width + (center.x + radius - x)] = color;

        if (error >= 0) {
            y += 1;
            error += 2 * (dy - dx);
        } else {
            error += 2 * dy;
        }
    }
}

void fill_hexagon(Image *image, isizex2 center, isize side, GifColorIndex color) {
    // Even side length => left/right sides take up an odd number of pixels => pixel perfect center.
    assert(side % 2 == 0);

    isize radius = (isize)roundf((f32)side * sqrtf(3.0F) / 2.0F);

    // Draw west and east sides and fill the space inbetween.
    for (isize y = center.y - side / 2; y <= center.y + side / 2; y += 1) {
        for (isize x = center.x - radius; x <= center.x + radius; x += 1) {
            image->pixels[y * image->width + x] = color;
        }
    }

    // Fill in the leftover south and north triangles using Bresenham.
    isize dx = radius;
    isize dy = side / 2;
    isize error = 2 * dy - dx;

    isize y = 0;
    for (isize x = 0; x <= radius; x += 1) {
        for (isize fill_x = center.x - radius + x; fill_x <= center.x + radius - x; fill_x += 1) {
            image->pixels[(center.y - side / 2 - y) * image->width + fill_x] = color;
            image->pixels[(center.y + side / 2 + y) * image->width + fill_x] = color;
        }

        if (error >= 0) {
            y += 1;
            error += 2 * (dy - dx);
        } else {
            error += 2 * dy;
        }
    }
}

typedef struct {
    bool *data;
    isize width;
    isize height;
    isize stride;
} Grid;

Grid sub_grid(Grid *grid, isize x, isize y, isize width, isize height) {
    return (Grid){
        .data = &grid->data[y * grid->stride + x],
        .width = width,
        .height = height,
        .stride = grid->stride,
    };
}

void grid_clear(Grid *grid) {
    for (isize grid_y = 0; grid_y < grid->height; grid_y += 1) {
        for (isize grid_x = 0; grid_x < grid->width; grid_x += 1) {
            grid->data[(grid->height - 1) * grid->stride + (grid->width - 1)] = false;
        }
    }
}

void grid_copy(Grid *source_grid, Grid *dest_grid) {
    assert(source_grid->width == dest_grid->width && source_grid->height == dest_grid->height);

    for (isize grid_y = 0; grid_y < source_grid->height; grid_y += 1) {
        memcpy(
            &dest_grid->data[grid_y * dest_grid->stride],
            &source_grid->data[grid_y * source_grid->stride],
            (size_t)(source_grid->width * sizeof(bool))
        );
    }
}

void grid_init_with_random(Grid *grid, PCG32 *rng) {
    for (isize grid_y = 0; grid_y < grid->height; grid_y += 1) {
        for (isize grid_x = 0; grid_x < grid->width; grid_x += 1) {
            u32 random = pcg32_random(rng);

            if (random < ((u32)1 << 28)) {
                grid->data[grid_y * grid->stride + grid_x] = true;
            } else {
                grid->data[grid_y * grid->stride + grid_x] = false;
            }
        }
    }
}

void grid_init_rule_110(Grid *grid) {
    grid_clear(grid);
    grid->data[(grid->height - 1) * grid->stride + (grid->width - 1)] = true;
}

void grid_step_rule_110(Grid *source_grid, Grid *dest_grid) {
    assert(source_grid->width == dest_grid->width && source_grid->height == dest_grid->height);
    assert(source_grid->height >= 2);

    isize grid_width = source_grid->width;
    isize grid_height = source_grid->height;

    for (isize grid_y = 0; grid_y < grid_height - 1; grid_y += 1) {
        memcpy(
            &dest_grid->data[grid_y * dest_grid->stride],
            &source_grid->data[(grid_y + 1) * source_grid->stride],
            (size_t)(grid_width * sizeof(bool))
        );
    }

    for (isize grid_x = 0; grid_x < grid_width; grid_x += 1) {
        bool left = dest_grid->data[
            (grid_height - 2) * dest_grid->stride +
            (grid_x - 1 + grid_width) % grid_width
        ];
        bool top = dest_grid->data[
            (grid_height - 2) * dest_grid->stride +
            grid_x
        ];
        bool right = dest_grid->data[
            (grid_height - 2) * dest_grid->stride +
            (grid_x + 1) % grid_width
        ];

        u32 pattern = (left ? (1 << 2) : 0) | (top ? (1 << 1) : 0) | (right ? 1 : 0);
        assert(0 <= pattern && pattern <= 7);

        bool transitions[] = {
            [7 /* 111 */] = 0,
            [6 /* 110 */] = 1,
            [5 /* 101 */] = 1,
            [4 /* 100 */] = 0,
            [3 /* 011 */] = 1,
            [2 /* 010 */] = 1,
            [1 /* 001 */] = 1,
            [0 /* 000 */] = 0,
        };

        dest_grid->data[(grid_height - 1) * dest_grid->stride + grid_x] = transitions[pattern];
    }
}

void grid_step_hex_automata(Grid *source_grid, Grid *dest_grid, bool wrap) {
    assert(source_grid->width == dest_grid->width && source_grid->height == dest_grid->height);

    isize grid_width = source_grid->width;
    isize grid_height = source_grid->height;

    for (isize grid_y = 0; grid_y < grid_height; grid_y += 1) {
        for (isize grid_x = 0; grid_x < grid_width; grid_x += 1) {
            isize dxs[] = {
                -1,                         // west
                1,                          // east
                grid_y % 2 == 0 ? -1 : 0,   // north-west
                grid_y % 2 == 0 ? 0 : 1,    // north-east
                grid_y % 2 == 0 ? -1 : 0,   // south-west
                grid_y % 2 == 0 ? 0 : 1,    // south-east
            };

            isize dys[] = {
                0,  // west
                0,  // east
                -1, // north-west
                -1, // north-east
                1,  // south-west
                1,  // south-east
            };

            isize alive_neighbors = 0;

            for (isize i = 0; i < 6; i += 1) {
                isize neighbor_x;
                isize neighbor_y;
                if (wrap) {
                    neighbor_x = (grid_x + dxs[i] + grid_width) % grid_width;
                    neighbor_y = (grid_y + dys[i] + grid_height) % grid_height;
                } else {
                    neighbor_x = grid_x + dxs[i];
                    neighbor_y = grid_y + dys[i];
                    if (
                        neighbor_x < 0 || neighbor_x >= grid_width ||
                        neighbor_y < 0 || neighbor_y >= grid_height
                    ) {
                        continue;
                    }
                }

                if (source_grid->data[neighbor_y * grid_width + neighbor_x]) {
                    alive_neighbors += 1;
                }
            }

            if (source_grid->data[grid_y * source_grid->stride + grid_x]) {
                if (alive_neighbors == 2) {
                    dest_grid->data[grid_y * dest_grid->stride + grid_x] = true;
                } else {
                    dest_grid->data[grid_y * dest_grid->stride + grid_x] = false;
                }
            } else {
                if (alive_neighbors == 2) {
                    dest_grid->data[grid_y * dest_grid->stride + grid_x] = true;
                } else {
                    dest_grid->data[grid_y * dest_grid->stride + grid_x] = false;
                }
            }
        }
    }
}

#define ARENA_CAPACITY (128 * 1024 * 1024)
typedef GifArena Arena;
#define arena_alloc gif_arena_alloc

typedef enum {
    AUTOHEX_MODE_FEED_110,
    AUTOHEX_MODE_RANDOM,
} AutohexMode;

int main(void) {
    PCG32 rng;
    pcg32_init(&rng, (u64)time(NULL));

    u8 *arena_memory = malloc(ARENA_CAPACITY);
    if (arena_memory == NULL) {
        return 1;
    }
    Arena arena = {arena_memory, arena_memory + ARENA_CAPACITY};

    AutohexMode mode = AUTOHEX_MODE_FEED_110;

    isize step_count = 250;

    isize grid_width = 80;
    isize grid_height = 80;
    isize hexagon_side = 8;
    isize hexagon_radius = (isize)roundf((f32)hexagon_side * sqrtf(3.0F) / 2.0F);

    Image image = {
        .width = grid_width * (2 * hexagon_radius) + (hexagon_radius + 1),
        .height = grid_height * (3 * hexagon_side / 2) + (hexagon_side / 2 + 1),
    };
    image.pixels = arena_alloc(&arena, image.width * image.height);

    isize color_count = srgb_palette_black_and_white(NULL);
    u8 *colors = arena_alloc(&arena, color_count * 3);
    srgb_palette_black_and_white(colors);

    FILE *output_file = fopen("autohex.gif", "wb");
    if (output_file == NULL) {
        return 1;
    }

    GifEncoder *encoder = gif_encoder_create(&arena);
    gif_frame_delay(encoder, 0.1F);

    GifOutputBuffer out_buffer;
    gif_out_buffer_create(8192, &out_buffer, &arena);

    Grid *grid_ping = arena_alloc(&arena, sizeof(Grid));
    grid_ping->data = arena_alloc(&arena, grid_width * grid_height * sizeof(bool));
    grid_ping->width = grid_width;
    grid_ping->height = grid_height;
    grid_ping->stride = grid_width;

    switch (mode) {
    case AUTOHEX_MODE_FEED_110: {
        grid_init_rule_110(grid_ping);
    } break;

    case AUTOHEX_MODE_RANDOM: {
        grid_init_with_random(grid_ping, &rng);
    } break;
    }

    Grid *grid_pong = arena_alloc(&arena, sizeof(Grid));
    grid_pong->data = arena_alloc(&arena, grid_width * grid_height * sizeof(bool));
    grid_pong->width = grid_width;
    grid_pong->height = grid_height;
    grid_pong->stride = grid_width;

    gif_encoder_start(encoder, image.width, image.height, colors, color_count, &out_buffer);
    if (gif_out_buffer_capacity_left(&out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
        fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
        gif_out_buffer_reset(&out_buffer);
    }

    for (isize step = 0; step < step_count; step += 1) {
        // Render the frame

        memset(image.pixels, 0, (size_t)(image.width * image.height));

        isize y = hexagon_side;
        for (isize grid_y = 0; grid_y < grid_height; grid_y += 1) {
            isize x;
            if (grid_y % 2 == 0) {
                x = hexagon_radius;
            } else {
                x = 2 * hexagon_radius;
            }

            for (isize grid_x = 0; grid_x < grid_width; grid_x += 1) {
                if (grid_ping->data[grid_y * grid_width + grid_x]) {
                    fill_hexagon(&image, (isizex2){x, y}, hexagon_side, 1);
                } else {
                    draw_hexagon(&image, (isizex2){x, y}, hexagon_side, 1);
                }

                x += 2 * hexagon_radius;
            }

            y += 3 * hexagon_side / 2;
        }

        printf("Frame %ld/%ld rendered\n", step + 1, step_count);

        // Encode the frame.

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

        printf("Frame %ld/%ld encoded\n", step + 1, step_count);

        // Advance the cellular automata state.

        switch (mode) {
        case AUTOHEX_MODE_FEED_110: {
            Grid bottom_ping = sub_grid(
                grid_ping,
                0, grid_height / 2,
                grid_width, grid_height - grid_height / 2
            );
            Grid bottom_pong = sub_grid(
                grid_pong,
                0, grid_height / 2,
                grid_width, grid_height - grid_height / 2
            );
            grid_step_rule_110(&bottom_ping, &bottom_pong);
            grid_copy(&bottom_pong, &bottom_ping);

            Grid top_ping = sub_grid(
                grid_ping,
                0, 0,
                grid_width, grid_height / 2 + 1
            );
            Grid top_pong = sub_grid(
                grid_pong,
                0, 0,
                grid_width, grid_height / 2 + 1
            );
            grid_step_hex_automata(&top_ping, &top_pong, false);
        } break;

        case AUTOHEX_MODE_RANDOM: {
            grid_step_hex_automata(grid_ping, grid_pong, true);
        } break;
        }

        Grid *swap = grid_ping;
        grid_ping = grid_pong;
        grid_pong = swap;
    }
    gif_encoder_finish(encoder, &out_buffer);

    // Write any leftover buffered data into the file.
    fwrite(out_buffer.data, 1, (size_t)out_buffer.encoded_size, output_file);
    fclose(output_file);

    return 0;
}
