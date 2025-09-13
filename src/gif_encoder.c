// GIF references:
// - https://www.w3.org/Graphics/GIF/spec-gif89a.txt
// - https://en.wikipedia.org/wiki/GIF
// - https://web.archive.org/web/20180620143135/https://www.vurdalakov.net/misc/gif/netscape-looping-application-extension

// TODO: Measure how much time (and memory) does each stage of converting an image to GIF take.
// TODO: Split hashset into multiple smaller ones when searching for unique colors in median-cut?
// TODO: Transdiff from ffmpeg: unchanged pixels are encoded as transparent on next frame.
// TODO: Color search: partition space, create list of possible closest colors for each partition.
// TODO: Color search: sort colors along 3 axis, first check colors which have closer projections.
// TODO: Try using K-means clustering instead of median-cut?
// TODO: Implement dithering? I don't know how it's done (though, random one should be easy).
// TODO: Pass the number of color components in arguments (1 for monochrome, 3 for RGB, 4 for RGBA).

#include "gif_encoder.h"

#include <assert.h> // assert
#include <stdlib.h> // qsort
#include <stddef.h> // size_t, NULL
#include <math.h>   // copysignf, powf, roundf, isnan, nanf, cbrtf, INFINITY
#include <string.h> // memmove, memcpy, memset, memcmp

#define sizeof(expr) (isize)sizeof(expr)
#define lengthof(string) (sizeof(string) - 1)
#define countof(expr) (sizeof(expr) / sizeof((expr)[0]))

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

static uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

static f64 f64_random(pcg32_random_t *rng) {
    // Multiply by 2^(-32).
    return (f64)pcg32_random_r(rng) * 0x1p-32;
}

#define ARENA_ALIGNMENT 16

typedef struct {
    u8 *begin;
    u8 *end;
} Arena;

// Just the "begin" pointer. The "end" pointer never changes after an allocation anyway.
typedef u8 *ArenaSnapshot;

static inline ArenaSnapshot arena_snapshot(Arena *arena) {
    return arena->begin;
}

static inline void arena_rewind(Arena *arena, ArenaSnapshot snapshot) {
    arena->begin = snapshot;
}

static inline isize arena_memory_left(Arena *arena) {
    return arena->end - arena->begin;
}

// Optional user-defined pre-allocation hook.
// Allows the user to define an out-of-memory behavior, as well as tracking the max memory usage.
//
// The user can also add more members to the arena struct if needed. Technically, I think this would
// violate strict aliasing, because then definitions of the arena struct would differ here in
// library code and in user code (despite both of them having the same starting sequence of fields).
//
// I believe there is an ugly way to avoid the UB (remove Arena struct definition from the library
// code and access the "begin" and "end" fields via memcpy), but I don't want to bother with it.

#ifndef GIF_LIB_REPORT_ALLOC
    #define GIF_LIB_REPORT_ALLOC(arena, size)
#else
    void GIF_LIB_REPORT_ALLOC(void *arena, isize size);
#endif

static void *arena_alloc(Arena *arena, isize size) {
    // Avoid reporting allocations on the "fake" arena which is used for calculating memory usage.
    const uptr UPTR_MAX = ~(uptr)0;
    if ((uptr)arena->end != UPTR_MAX) {
        GIF_LIB_REPORT_ALLOC(arena, size);
    }

    if (size == 0) {
        return NULL;
    }

    isize padding = (~(uptr)arena->begin + 1) & (ARENA_ALIGNMENT - 1);
    isize memory_left = arena->end - arena->begin - padding;
    if (memory_left < size) {
        return NULL;
    }

    void *ptr = arena->begin + padding;
    arena->begin += padding + size;
    return ptr;
}

static void *arena_realloc(Arena *arena, void *old_ptr, isize old_size, isize new_size) {
    assert(old_ptr != NULL);
    assert(((uptr)old_ptr & (ARENA_ALIGNMENT - 1)) == 0);

    if (new_size == 0) {
        return NULL;
    }
    if (old_size >= new_size) {
        return old_ptr;
    }

    if ((u8 *)old_ptr + old_size == arena->begin) {
        arena->begin = old_ptr;
        return arena_alloc(arena, new_size);
    } else {
        void *new_ptr = arena_alloc(arena, new_size);
        if (new_ptr == NULL) {
            return NULL;
        }

        memcpy(new_ptr, old_ptr, (size_t)old_size);
        return new_ptr;
    }
}

#define COMPONENTS_PER_COLOR 3

u8 *srgb_palette_black_and_white(isize *color_count, void *arena) {
    u8 *colors = arena_alloc(arena, 2 * COMPONENTS_PER_COLOR * sizeof(u8));
    if (colors == NULL) {
        *color_count = 0;
        return NULL;
    }

    u8 *color_iter = colors;

    color_iter[0] = 0x00;
    color_iter[1] = 0x00;
    color_iter[2] = 0x00;
    color_iter += COMPONENTS_PER_COLOR;

    color_iter[0] = 0xff;
    color_iter[1] = 0xff;
    color_iter[2] = 0xff;
    color_iter += COMPONENTS_PER_COLOR;

    *color_count = 2;
    return colors;
}

u8 *srgb_palette_monochrome(isize *color_count, void *arena) {
    u8 *colors = arena_alloc(arena, 256 * COMPONENTS_PER_COLOR * sizeof(u8));
    if (colors == NULL) {
        *color_count = 0;
        return NULL;
    }

    u8 *color_iter = colors;

    for (int shade = 0x00; shade <= 0xff; shade += 1) {
        color_iter[0] = (u8)shade;
        color_iter[1] = (u8)shade;
        color_iter[2] = (u8)shade;

        color_iter += COMPONENTS_PER_COLOR;
    }

    *color_count = 256;
    return colors;
}

u8 *srgb_palette_web_safe(isize *color_count, void *arena) {
    // 216 colors, because 6 * 6 * 6 == 216.
    u8 *colors = arena_alloc(arena, 216 * COMPONENTS_PER_COLOR * sizeof(u8));
    if (colors == NULL) {
        *color_count = 0;
        return NULL;
    }

    u8 *color_iter = colors;

    for (int red = 0x00; red <= 0xff; red += 0x33) {
        for (int green = 0x00; green <= 0xff; green += 0x33) {
            for (int blue = 0x00; blue <= 0xff; blue += 0x33) {
                color_iter[0] = (u8)red;
                color_iter[1] = (u8)green;
                color_iter[2] = (u8)blue;

                color_iter += COMPONENTS_PER_COLOR;
            }
        }
    }

    *color_count = 216;
    return colors;
}

f32 *srgb_to_float(u8 const *srgb_colors, isize color_count, void *arena) {
    f32 *float_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (float_colors == NULL) {
        return NULL;
    }

    f32 *float_color_iter = float_colors;

    u8 const *srgb_color_iter = srgb_colors;

    while (srgb_color_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        float_color_iter[0] = (f32)srgb_color_iter[0] / 255.0F;
        float_color_iter[1] = (f32)srgb_color_iter[1] / 255.0F;
        float_color_iter[2] = (f32)srgb_color_iter[2] / 255.0F;

        srgb_color_iter += COMPONENTS_PER_COLOR;
        float_color_iter += COMPONENTS_PER_COLOR;
    }

    return float_colors;
}

u8 *float_to_srgb(f32 const *float_colors, isize color_count, void *arena) {
    u8 *srgb_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (srgb_colors == NULL) {
        return NULL;
    }

    u8 *srgb_color_iter = srgb_colors;

    f32 const *float_color_iter = float_colors;

    while (float_color_iter != float_colors + color_count * COMPONENTS_PER_COLOR) {
        srgb_color_iter[0] = (u8)(float_color_iter[0] * 255.0F);
        srgb_color_iter[1] = (u8)(float_color_iter[1] * 255.0F);
        srgb_color_iter[2] = (u8)(float_color_iter[2] * 255.0F);

        float_color_iter += COMPONENTS_PER_COLOR;
        srgb_color_iter += COMPONENTS_PER_COLOR;
    }

    return srgb_colors;
}

// https://en.wikipedia.org/wiki/SRGB#Transfer_function_("gamma")
static inline f32 srgb_component_to_linear(f32 c) {
    return c <= 0.04045F ? c / 12.92F : powf((c + 0.055F) / 1.055F, 2.4F);
}

f32 *srgb_to_linear(u8 const *srgb_colors, isize color_count, void *arena) {
    f32 *linear_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (linear_colors == NULL) {
        return NULL;
    }

    f32 *linear_color_iter = linear_colors;

    u8 const *srgb_color_iter = srgb_colors;

    while (srgb_color_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        linear_color_iter[0] = srgb_component_to_linear((f32)srgb_color_iter[0] / 255.0F);
        linear_color_iter[1] = srgb_component_to_linear((f32)srgb_color_iter[1] / 255.0F);
        linear_color_iter[2] = srgb_component_to_linear((f32)srgb_color_iter[2] / 255.0F);

        srgb_color_iter += COMPONENTS_PER_COLOR;
        linear_color_iter += COMPONENTS_PER_COLOR;
    }

    return linear_colors;
}

// https://en.wikipedia.org/wiki/SRGB#Transfer_function_("gamma")
static inline f32 linear_component_to_srgb(f32 c) {
    return c <= 0.0031308F ? 12.92F * c : 1.055F * powf(c, 1.0F / 2.4F) - 0.055F;
}

u8 *linear_to_srgb(f32 const *linear_colors, isize color_count, void *arena) {
    u8 *srgb_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (srgb_colors == NULL) {
        return NULL;
    }

    u8 *srgb_color_iter = srgb_colors;

    f32 const *linear_color_iter = linear_colors;

    while (srgb_color_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        srgb_color_iter[0] = (u8)(linear_component_to_srgb(linear_color_iter[0]) * 255.0F);
        srgb_color_iter[1] = (u8)(linear_component_to_srgb(linear_color_iter[1]) * 255.0F);
        srgb_color_iter[2] = (u8)(linear_component_to_srgb(linear_color_iter[2]) * 255.0F);

        srgb_color_iter += COMPONENTS_PER_COLOR;
        linear_color_iter += COMPONENTS_PER_COLOR;
    }

    return srgb_colors;
}

typedef struct {
    f32 x;
    f32 y;
    f32 z;
} f32x3;

f32 *srgb_to_lab(u8 const *srgb_colors, isize color_count, void *arena) {
    f32 *lab_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (lab_colors == NULL) {
        return NULL;
    }

    f32 *lab_color_iter = lab_colors;

    u8 const *srgb_color_iter = srgb_colors;

    while (srgb_color_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        f32x3 linear_srgb = {
            srgb_component_to_linear((f32)srgb_color_iter[0] / 255.0F) * 100.0F,
            srgb_component_to_linear((f32)srgb_color_iter[1] / 255.0F) * 100.0F,
            srgb_component_to_linear((f32)srgb_color_iter[2] / 255.0F) * 100.0F,
        };

        // sRGB -> XYZ (D65)
        // https://en.wikipedia.org/wiki/SRGB#Primaries
        f32x3 xyz = {
            0.4124F * linear_srgb.x + 0.3576F * linear_srgb.y + 0.1805F * linear_srgb.z,
            0.2126F * linear_srgb.x + 0.7152F * linear_srgb.y + 0.0722F * linear_srgb.z,
            0.0193F * linear_srgb.x + 0.1192F * linear_srgb.y + 0.9505F * linear_srgb.z,
        };

        // Standard Illuminant D65
        f32 Xr = 95.0489F;
        f32 Yr = 100.0F;
        f32 Zr = 108.8840F;

        // https://en.wikipedia.org/wiki/CIELAB_color_space
        //
        // The lightness value, L*, defines black at 0 and white at 100.
        //
        // The a* axis is relative to the green–red opponent colors, with negative values toward
        // green and positive values toward red. The b* axis represents the blue–yellow opponents,
        // with negative numbers toward blue and positive toward yellow.
        //
        // The a* and b* axes are unbounded and depending on the reference white they can easily
        // exceed ±150 to cover the human gamut. Nevertheless, software implementations often clamp
        // these values for practical reasons. For instance, if integer math is being used it is
        // common to clamp a* and b* in the range of −128 to 127.

        // https://en.wikipedia.org/wiki/CIELAB_color_space#From_CIE_XYZ_to_CIELAB

        f32 xr = xyz.x / Xr;
        f32 fx = xr > 216.0F / 24389.0F ? cbrtf(xr) : (24389.0F * xr / 27.0F + 16.0F) / 116.0F;

        f32 yr = xyz.y / Yr;
        f32 fy = yr > 216.0F / 24389.0F ? cbrtf(yr) : (24389.0F * yr / 27.0F + 16.0F) / 116.0F;

        f32 zr = xyz.z / Zr;
        f32 fz = zr > 216.0F / 24389.0F ? cbrtf(zr) : (24389.0F * zr / 27.0F + 16.0F) / 116.0F;

        int const L = 0, a = 1, b = 2;
        lab_color_iter[L] = 116.0F * fy - 16.0F;
        lab_color_iter[a] = 500.0F * (fx - fy);
        lab_color_iter[b] = 200.0F * (fy - fz);

        srgb_color_iter += COMPONENTS_PER_COLOR;
        lab_color_iter += COMPONENTS_PER_COLOR;
    }

    return lab_colors;
}

u8 *lab_to_srgb(f32 const *lab_colors, isize color_count, void *arena) {
    u8 *srgb_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (srgb_colors == NULL) {
        return NULL;
    }

    u8 *srgb_color_iter = srgb_colors;

    f32 const *lab_color_iter = lab_colors;

    while (lab_color_iter != lab_colors + color_count * COMPONENTS_PER_COLOR) {
        // Standard Illuminant D65
        f32 Xr = 95.0489F;
        f32 Yr = 100.0F;
        f32 Zr = 108.8840F;

        // https://en.wikipedia.org/wiki/CIELAB_color_space#From_CIELAB_to_CIEXYZ

        int const L = 0, a = 1, b = 2;

        f32 fy = (lab_color_iter[L] + 16.0F) / 116.0F;
        f32 fx = lab_color_iter[a] / 500.0F + fy;
        f32 fz = fy - lab_color_iter[b] / 200.0F;

        f32 xr = fx * fx * fx > 216.0F / 24389.0F
            ? fx * fx * fx
            : (116.0F * fx - 16.0F) * 27.0F / 24389.0F;

        f32 yr = lab_color_iter[L] > 216.0F / 27.0F
            ? fy * fy * fy
            : lab_color_iter[L] * 27.0F / 24389.0F;

        f32 zr = fz * fz * fz > 216.0F / 24389.0F
            ? fz * fz * fz
            : (116.0F * fz - 16.0F) * 27.0F / 24389.0F;

        f32x3 xyz = {xr * Xr, yr * Yr, zr * Zr};

        // XYZ -> sRGB (D65)
        // https://en.wikipedia.org/wiki/SRGB#Primaries
        f32x3 linear_srgb = {
            (+3.2406255F) * xyz.x + (-1.5372080F) * xyz.y + (-0.4986286F) * xyz.z,
            (-0.9689307F) * xyz.x + (+1.8757561F) * xyz.y + (+0.0415175F) * xyz.z,
            (+0.0557101F) * xyz.x + (-0.2040211F) * xyz.y + (+1.0569959F) * xyz.z,
        };

        srgb_color_iter[0] = (u8)(linear_component_to_srgb(linear_srgb.x / 100.0F) * 255.0F);
        srgb_color_iter[1] = (u8)(linear_component_to_srgb(linear_srgb.y / 100.0F) * 255.0F);
        srgb_color_iter[2] = (u8)(linear_component_to_srgb(linear_srgb.z / 100.0F) * 255.0F);

        lab_color_iter += COMPONENTS_PER_COLOR;
        srgb_color_iter += COMPONENTS_PER_COLOR;
    }

    return srgb_colors;
}

// https://en.wikipedia.org/wiki/Oklab_color_space#Conversion_from_sRGB
// https://bottosson.github.io/posts/oklab/#converting-from-linear-srgb-to-oklab
f32 *srgb_to_oklab(u8 const *srgb_colors, isize color_count, void *arena) {
    f32 *oklab_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (oklab_colors == NULL) {
        return NULL;
    }

    f32 *oklab_iter = oklab_colors;

    u8 const *srgb_iter = srgb_colors;

    while (srgb_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        f32x3 linear = {
            srgb_component_to_linear((f32)srgb_iter[0] / 255.0F),
            srgb_component_to_linear((f32)srgb_iter[1] / 255.0F),
            srgb_component_to_linear((f32)srgb_iter[2] / 255.0F),
        };

        // From Wikipedia:
        // > The (l,m,s) space used here is not the same as the LMS color space, but rather an
        // > arbitrary space that was found numerically to best fit the color appearance data.
        f32x3 lms = {
            0.4122214708f * linear.x + 0.5363325363f * linear.y + 0.0514459929f * linear.z,
            0.2119034982f * linear.x + 0.6806995451f * linear.y + 0.1073969566f * linear.z,
            0.0883024619f * linear.x + 0.2817188376f * linear.y + 0.6299787005f * linear.z,
        };
        lms = (f32x3){cbrtf(lms.x), cbrtf(lms.y), cbrtf(lms.z)};

        int const L = 0, a = 1, b = 2;
        oklab_iter[L] = 0.2104542553f * lms.x + 0.7936177850f * lms.y - 0.0040720468f * lms.z;
        oklab_iter[a] = 1.9779984951f * lms.x - 2.4285922050f * lms.y + 0.4505937099f * lms.z;
        oklab_iter[b] = 0.0259040371f * lms.x + 0.7827717662f * lms.y - 0.8086757660f * lms.z;

        srgb_iter += COMPONENTS_PER_COLOR;
        oklab_iter += COMPONENTS_PER_COLOR;
    }

    return oklab_colors;
}

// See srgb_to_oklab for the references.
u8 *oklab_to_srgb(f32 const *oklab_colors, isize color_count, void *arena) {
    u8 *srgb_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (srgb_colors == NULL) {
        return NULL;
    }

    u8 *srgb_iter = srgb_colors;

    f32 const *oklab_iter = oklab_colors;

    while (oklab_iter != oklab_colors + color_count * COMPONENTS_PER_COLOR) {
        int const L = 0, a = 1, b = 2;
        f32x3 lms = {
            oklab_iter[L] + 0.3963377774f * oklab_iter[a] + 0.2158037573f * oklab_iter[b],
            oklab_iter[L] - 0.1055613458f * oklab_iter[a] - 0.0638541728f * oklab_iter[b],
            oklab_iter[L] - 0.0894841775f * oklab_iter[a] - 1.2914855480f * oklab_iter[b],
        };
        lms = (f32x3){
            lms.x * lms.x * lms.x,
            lms.y * lms.y * lms.y,
            lms.z * lms.z * lms.z,
        };

        f32x3 linear = {
            +4.0767416621f * lms.x - 3.3077115913f * lms.y + 0.2309699292f * lms.z,
            -1.2684380046f * lms.x + 2.6097574011f * lms.y - 0.3413193965f * lms.z,
            -0.0041960863f * lms.x - 0.7034186147f * lms.y + 1.7076147010f * lms.z,
        };

        srgb_iter[0] = (u8)(linear_component_to_srgb(linear.x) * 255.0F);
        srgb_iter[1] = (u8)(linear_component_to_srgb(linear.y) * 255.0F);
        srgb_iter[2] = (u8)(linear_component_to_srgb(linear.z) * 255.0F);

        oklab_iter += COMPONENTS_PER_COLOR;
        srgb_iter += COMPONENTS_PER_COLOR;
    }

    return srgb_colors;
}

// Based on djb2 hashing function:
// http://www.cse.yorku.ca/~oz/hash.html
static u32 f32x3_hash(f32x3 value) {
    u32 hash = 5381;

    u8 const *value_iter = (u8 const *)&value;
    u8 const *value_end = (u8 const *)&value + sizeof(f32x3);
    while (value_iter != value_end) {
        hash = ((hash << 5) + hash) + *value_iter;
        value_iter += 1;
    }

    return hash;
}

static inline bool f32x3_equals(f32x3 left, f32x3 right) {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

// First color component set to NaN indicates an empty bucket.
typedef f32x3 ColorSetBucket;

static inline bool color_set_bucket_empty(ColorSetBucket bucket) {
    return isnan(bucket.x);
}

typedef struct {
    ColorSetBucket *buckets;
    isize bucket_count;
    isize color_count;
} ColorSet;

static void color_set_insert(ColorSet *colors, f32x3 new_color) {
    assert(colors->color_count < colors->bucket_count);

    isize bucket_index = f32x3_hash(new_color) % colors->bucket_count;
    while (
        !color_set_bucket_empty(colors->buckets[bucket_index]) &&
        !f32x3_equals(colors->buckets[bucket_index], new_color)
    ) {
        bucket_index = (bucket_index + 1) % colors->bucket_count;
    }

    if (color_set_bucket_empty(colors->buckets[bucket_index])) {
        colors->buckets[bucket_index] = new_color;
        colors->color_count += 1;
    }
}

static inline isize isize_min(isize left, isize right) {
    return left < right ? left : right;
}

static inline isize isize_max(isize left, isize right) {
    return left > right ? left : right;
}

static inline isize isize_abs(isize value) {
    if (value >= 0) {
        return value;
    } else {
        return -value;
    }
}

static int f32x3_compare_x_only(void const *left_ptr, void const *right_ptr) {
    f32x3 left = *(f32x3 *)left_ptr;
    f32x3 right = *(f32x3 *)right_ptr;

    return (int)copysignf(1.0F, left.x - right.x);
}

static int f32x3_compare_y_only(void const *left_ptr, void const *right_ptr) {
    f32x3 left = *(f32x3 *)left_ptr;
    f32x3 right = *(f32x3 *)right_ptr;

    return (int)copysignf(1.0F, left.y - right.y);
}

static int f32x3_compare_z_only(void const *left_ptr, void const *right_ptr) {
    f32x3 left = *(f32x3 *)left_ptr;
    f32x3 right = *(f32x3 *)right_ptr;

    return (int)copysignf(1.0F, left.z - right.z);
}

f32 *palette_by_median_cut(
    f32 const *pixels, isize pixel_count,
    isize target_color_count,
    isize *color_count,
    void *arena
) {
    assert(target_color_count <= GIF_MAX_COLORS);

    f32 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (palette == NULL) {
        *color_count = 0;
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    // TODO: I don't like this hashset capacity estimate, the test image only uses around 2% of the
    // allocated buckets. And also keep in mind that you will have to iterate over the entire list
    // of buckets to extract the unique colors. Maybe rehashing here is not a bad idea after all?

    ColorSet color_set;
    color_set.color_count = 0;

    // Worst case scenario: all pixels in the input image are different.
    color_set.bucket_count = isize_min(pixel_count, 256 * 256 * 256) * 3 / 2;
    color_set.buckets = arena_alloc(arena, color_set.bucket_count * sizeof(ColorSetBucket));
    if (color_set.buckets == NULL) {
        *color_count = 0;
        return NULL;
    }

    // If the implementation does not support quiet NaNs, these functions return zero.
    assert(nanf("") != 0.0F);
    for (isize i = 0; i < color_set.bucket_count; i += 1) {
        color_set.buckets[i].x = nanf("");
        color_set.buckets[i].y = nanf("");
        color_set.buckets[i].z = nanf("");
    }

    f32 const *pixel_iter = pixels;
    f32 const *pixels_end = pixels + COMPONENTS_PER_COLOR * pixel_count;
    while (pixel_iter < pixels_end) {
        color_set_insert(&color_set, (f32x3){pixel_iter[0], pixel_iter[1], pixel_iter[2]});
        pixel_iter += COMPONENTS_PER_COLOR;
    }

    f32x3 *colors = arena_alloc(arena, color_set.color_count * sizeof(f32x3));
    if (colors == NULL) {
        *color_count = 0;
        return NULL;
    }

    f32x3 *color_iter = colors;
    f32x3 *colors_end = colors + color_set.color_count;

    ColorSetBucket *bucket_iter = color_set.buckets;
    while (bucket_iter < color_set.buckets + color_set.bucket_count) {
        if (!color_set_bucket_empty(*bucket_iter)) {
            assert(color_iter != colors_end);

            *color_iter = *bucket_iter;
            color_iter += 1;
        }

        bucket_iter += 1;
    }
    assert(color_iter == colors_end);

    typedef struct {
        f32x3 *begin;
        f32x3 *end;
    } Segment;

    // Ring buffer is used for the queue.
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
    queue.segments[queue.first_index] = (Segment){colors, colors_end};

    while (queue.count < target_color_count) {
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

        f32 min_x = INFINITY, min_y = INFINITY, min_z = INFINITY;
        f32 max_x = -INFINITY, max_y = -INFINITY, max_z = -INFINITY;
        {
            f32x3 *colors_iter = current_segment.begin;
            while (colors_iter != current_segment.end) {
                min_x = min_x < colors_iter->x ? min_x : colors_iter->x;
                max_x = max_x > colors_iter->x ? max_x : colors_iter->x;

                min_y = min_y < colors_iter->y ? min_y : colors_iter->y;
                max_y = max_y > colors_iter->y ? max_y : colors_iter->y;

                min_z = min_z < colors_iter->z ? min_z : colors_iter->z;
                max_z = max_z > colors_iter->z ? max_z : colors_iter->z;

                colors_iter += 1;
            }
        }

        int(*comparator)(void const *, void const *);
        if (max_x - min_x >= max_y - min_y && max_x - min_x >= max_z - min_z) {
            comparator = f32x3_compare_x_only;
        } else if (max_y - min_y >= max_x - min_x && max_y - min_y >= max_z - min_z) {
            comparator = f32x3_compare_y_only;
        } else {
            comparator = f32x3_compare_z_only;
        }
        qsort(
            current_segment.begin,
            (size_t)(current_segment.end - current_segment.begin),
            sizeof(f32x3),
            comparator
        );

        queue.segments[(queue.last_index + 1) % countof(queue.segments)] = (Segment){
            .begin = current_segment.begin,
            .end = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
        };
        queue.segments[(queue.last_index + 2) % countof(queue.segments)] = (Segment){
            .begin = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
            .end = current_segment.end,
        };
        queue.last_index = (queue.last_index + 2) % countof(queue.segments);
        queue.count += 2;
    }

    f32 *quantized_colors_iter = palette;
    isize segment_index = queue.first_index;

    while (true) {
        f64 x_sum = 0.0;
        f64 y_sum = 0.0;
        f64 z_sum = 0.0;

        f32x3 *segment_iter = queue.segments[segment_index].begin;
        while (segment_iter < queue.segments[segment_index].end) {
            x_sum += segment_iter->x;
            y_sum += segment_iter->y;
            z_sum += segment_iter->z;

            segment_iter += 1;
        }

        isize color_count = queue.segments[segment_index].end - queue.segments[segment_index].begin;
        quantized_colors_iter[0] = (f32)(x_sum / (f64)color_count);
        quantized_colors_iter[1] = (f32)(y_sum / (f64)color_count);
        quantized_colors_iter[2] = (f32)(z_sum / (f64)color_count);

        if (segment_index == queue.last_index) {
            break;
        }
        quantized_colors_iter += COMPONENTS_PER_COLOR;
        segment_index = (segment_index + 1) % countof(queue.segments);
    }

    arena_rewind(arena, snapshot);
    *color_count = queue.count;
    return palette;
}

f32 *palette_by_k_means(
    f32 const *pixels, isize pixel_count,
    isize target_color_count,
    isize *colors_generated,
    void *arena_void
) {
    isize const MAX_ITERATIONS = 300;

    Arena *arena = arena_void;

    f32 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (palette == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    // Free space left from deallocated hash sets.
    Arena free_space = {.begin = arena->begin, .end = arena->begin};

    ColorSet color_set = {0};
    color_set.bucket_count = 1024;
    color_set.buckets = arena_alloc(arena, 1024 * sizeof(ColorSetBucket));
    if (color_set.buckets == NULL) {
        *colors_generated = 0;
        return NULL;
    }
    for (isize i = 0; i < color_set.bucket_count; i += 1) {
        color_set.buckets[i].x = nanf("");
    }

    f32 const *pixel_iter = pixels;
    f32 const *pixels_end = pixels + COMPONENTS_PER_COLOR * pixel_count;
    while (pixel_iter < pixels_end) {
        if ((f32)color_set.color_count / color_set.bucket_count > 0.75F) {
            // The growth factor is set to 1.5 so that the free_space arena eventually grows large
            // enough to fit a new rehashed hashset.
            isize new_bucket_count = color_set.bucket_count * 3 / 2;

            ColorSetBucket *new_buckets = NULL;
            if (arena_memory_left(&free_space) >= new_bucket_count * sizeof(ColorSetBucket)) {
                new_buckets = (ColorSetBucket *)free_space.begin;
                arena->begin = free_space.end;
                free_space.end = free_space.begin;
            } else {
                new_buckets = arena_alloc(arena, new_bucket_count * sizeof(ColorSetBucket));
                free_space.end += color_set.bucket_count * sizeof(ColorSetBucket);
            }
            if (new_buckets == NULL) {
                *colors_generated = 0;
                return NULL;
            }

            for (isize i = 0; i < new_bucket_count; i += 1) {
                new_buckets[i].x = nanf("");
            }

            ColorSet new_color_set = {
                .buckets = new_buckets,
                .bucket_count = new_bucket_count,
                .color_count = 0,
            };
            for (isize i = 0; i < color_set.bucket_count; i += 1) {
                if (!color_set_bucket_empty(color_set.buckets[i])) {
                    color_set_insert(&new_color_set, color_set.buckets[i]);
                }
            }
            assert(color_set.color_count == new_color_set.color_count);

            color_set = new_color_set;
        }

        color_set_insert(&color_set, (f32x3){pixel_iter[0], pixel_iter[1], pixel_iter[2]});
        pixel_iter += COMPONENTS_PER_COLOR;
    }

    arena_rewind(arena, snapshot);
    f32x3 *colors = arena_alloc(arena, color_set.color_count * sizeof(f32x3));
    isize color_count = color_set.color_count;

    // We already have all of the colors in memory, so there is definitely enough space.
    assert(colors != NULL);

    {
        f32x3 *color_iter = colors;
        ColorSetBucket *bucket_iter = color_set.buckets;
        while (bucket_iter < color_set.buckets + color_set.bucket_count) {
            if (!color_set_bucket_empty(*bucket_iter)) {
                assert(color_iter != colors + color_count);

                *color_iter = *bucket_iter;
                color_iter += 1;
            }

            bucket_iter += 1;
        }
        assert(color_iter == colors + color_count);
    }

    // In total the whole image contained less colors than we asked for.
    if (color_count <= target_color_count) {
        memcpy(palette, colors, color_count * sizeof(f32x3));
        *colors_generated = color_count;
        arena_rewind(arena, snapshot);
        return palette;
    }

    // Pick the initial set of centroids randomly from the image colors.
    isize centroid_count = target_color_count;
    f32x3 *centroids = arena_alloc(arena, centroid_count * sizeof(f32x3));
    f32x3 *new_centroids = arena_alloc(arena, centroid_count * sizeof(f32x3));
    if (centroids == NULL || new_centroids == NULL) {
        *colors_generated = 0;
        return NULL;
    }
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (u64)memcpy, (u64)memmove);
    {
        f32x3 *centroid_iter = centroids;
        f32x3 *color_iter = colors, *colors_end = colors + color_count;
        for (isize i = 0; i < centroid_count; i += 1) {
            isize next_color_index = f64_random(&rng) * (colors_end - color_iter);
            assert(next_color_index < colors_end - color_iter);

            *centroid_iter = color_iter[next_color_index];
            centroid_iter += 1;

            f32x3 swap = color_iter[next_color_index];
            color_iter[next_color_index] = color_iter[0];
            color_iter[0] = swap;
            color_iter += 1;
        }
    }

    isize *colors_to_clusters = arena_alloc(arena, color_count * sizeof(isize));

    isize cluster_count = target_color_count;
    isize *cluster_sizes = arena_alloc(arena, cluster_count * sizeof(isize));

    if (colors_to_clusters == NULL || cluster_sizes == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    // Assign every single color to the 0th cluster initially.
    memset(colors_to_clusters, 0, color_count * sizeof(isize));
    memset(cluster_sizes, 0, cluster_count * sizeof(isize));
    cluster_sizes[0] = color_count;

    // Arrays used for averaging colors within the clusters to get the new centroids.
    f64 *x_sum = arena_alloc(arena, cluster_count * sizeof(f64));
    f64 *y_sum = arena_alloc(arena, cluster_count * sizeof(f64));
    f64 *z_sum = arena_alloc(arena, cluster_count * sizeof(f64));
    if (x_sum == NULL || y_sum == NULL || z_sum == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    for (isize iteration = 0; iteration < MAX_ITERATIONS; iteration += 1) {
        isize colors_which_changed_clusters = 0;

        for (isize i = 0; i < color_count; i += 1) {
            f32 min_distance =
                (colors[i].x - centroids[0].x) * (colors[i].x - centroids[0].x) +
                (colors[i].y - centroids[0].y) * (colors[i].y - centroids[0].y) +
                (colors[i].z - centroids[0].z) * (colors[i].z - centroids[0].z);
            isize closest_cluster_index = 0;

            for (isize j = 1; j < centroid_count; j += 1) {
                f32 distance =
                    (colors[i].x - centroids[j].x) * (colors[i].x - centroids[j].x) +
                    (colors[i].y - centroids[j].y) * (colors[i].y - centroids[j].y) +
                    (colors[i].z - centroids[j].z) * (colors[i].z - centroids[j].z);

                if (distance < min_distance) {
                    min_distance = distance;
                    closest_cluster_index = j;
                }
            }

            if (colors_to_clusters[i] != closest_cluster_index) {
                isize old_cluster = colors_to_clusters[i];
                cluster_sizes[old_cluster] -= 1;

                isize new_cluster = closest_cluster_index;
                cluster_sizes[new_cluster] += 1;

                colors_which_changed_clusters += 1;
                colors_to_clusters[i] = new_cluster;
            }
        }

        // Quit early once at least 99% of the colors stopped changing their clusters.
        if (colors_which_changed_clusters <= color_count / 100) {
            break;
        }

        memset(x_sum, 0, centroid_count * sizeof(f64));
        memset(y_sum, 0, centroid_count * sizeof(f64));
        memset(z_sum, 0, centroid_count * sizeof(f64));
        for (isize i = 0; i < color_count; i += 1) {
            x_sum[colors_to_clusters[i]] += colors[i].x;
            y_sum[colors_to_clusters[i]] += colors[i].y;
            z_sum[colors_to_clusters[i]] += colors[i].z;
        }
        for (isize i = 0; i < centroid_count; i += 1) {
            new_centroids[i].x = x_sum[i] / cluster_sizes[i];
            new_centroids[i].y = y_sum[i] / cluster_sizes[i];
            new_centroids[i].z = z_sum[i] / cluster_sizes[i];
        }

        f32x3 *swap = centroids;
        centroids = new_centroids;
        new_centroids = swap;
    }

    memcpy(palette, centroids, centroid_count * sizeof(f32x3));
    *colors_generated = centroid_count;
    arena_rewind(arena, snapshot);
    return palette;
}

typedef struct {
    f32x3 color;
    GifColorIndex color_index;
} IndexedColor;

typedef int (*IndexedColorComparator)(void const *, void const *);

static int indexed_color_compare_by_x(void const *left, void const *right) {
    f32x3 left_color = ((IndexedColor *)left)->color;
    f32x3 right_color = ((IndexedColor *)right)->color;

    return (int)copysignf(1.0F, left_color.x - right_color.x);
}

static int indexed_color_compare_by_y(void const *left, void const *right) {
    f32x3 left_color = ((IndexedColor *)left)->color;
    f32x3 right_color = ((IndexedColor *)right)->color;

    return (int)copysignf(1.0F, left_color.y - right_color.y);
}

static int indexed_color_compare_by_z(void const *left, void const *right) {
    f32x3 left_color = ((IndexedColor *)left)->color;
    f32x3 right_color = ((IndexedColor *)right)->color;

    return (int)copysignf(1.0F, left_color.z - right_color.z);
}

typedef struct ColorTreeNode ColorTreeNode;

// K-d tree for closest color lookups.
struct ColorTreeNode {
    f32 median;
    ColorTreeNode *left;
    ColorTreeNode *right;
    IndexedColor value;
};

static ColorTreeNode *color_tree_construct(
    IndexedColor *colors,
    isize color_count,
    isize depth,
    Arena *arena
) {
    assert(color_count != 0);

    ColorTreeNode *node = arena_alloc(arena, sizeof(ColorTreeNode));

    if (color_count > 1) {
        isize color_component_index = depth % 3;

        IndexedColorComparator comparator = (IndexedColorComparator[]){
            indexed_color_compare_by_x,
            indexed_color_compare_by_y,
            indexed_color_compare_by_z,
        }[color_component_index];
        qsort(colors, (size_t)color_count, sizeof(IndexedColor), comparator);

        isize middle = color_count / 2;

        IndexedColor median = colors[middle];
        node->median = (f32[]){
            median.color.x,
            median.color.y,
            median.color.z,
        }[color_component_index];

        // The median value goes into the right subtree.
        node->left = color_tree_construct(colors, middle, depth + 1, arena);
        node->right = color_tree_construct(colors + middle, color_count - middle, depth + 1, arena);
    } else {
        node->left = NULL;
        node->right = NULL;
        node->value = colors[0];
    }

    return node;
}

// distance_to_nearest is actually squared distance.
static void color_tree_get_nearest(
    ColorTreeNode *node,
    f32x3 search_color,
    isize depth,
    f32 *distance_to_nearest,
    GifColorIndex *nearest_color_index
) {
    if (node->left != NULL) {
        f32 search_color_component = (f32[]){
            search_color.x,
            search_color.y,
            search_color.z,
        }[depth % 3];

        // This one is checked always:
        ColorTreeNode *subtree;
        // This one is checked only if needed:
        ColorTreeNode *other_subtree;
        if (search_color_component < node->median) {
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
            distance_to_nearest,
            nearest_color_index
        );

        f32 other_subtree_min_distance =
            (node->median - search_color_component) * (node->median - search_color_component);
        if (other_subtree_min_distance < *distance_to_nearest) {
            color_tree_get_nearest(
                other_subtree,
                search_color,
                depth + 1,
                distance_to_nearest,
                nearest_color_index
            );
        }
    } else {
        // Leaf node
        assert(node->left == NULL && node->right == NULL);

        f32 current_distance =
            (node->value.color.x - search_color.x) * (node->value.color.x - search_color.x) +
            (node->value.color.y - search_color.y) * (node->value.color.y - search_color.y) +
            (node->value.color.z - search_color.z) * (node->value.color.z - search_color.z);

        if (current_distance < *distance_to_nearest) {
            *distance_to_nearest = current_distance;
            *nearest_color_index = node->value.color_index;
        }
    }
}

GifColorIndex *image_quantize_for_gif(
    f32 const *pixels, isize pixel_count,
    f32 const *colors, isize color_count,
    void *arena
) {
    assert(color_count <= (isize)GIF_COLOR_INDEX_MAX + 1);

    GifColorIndex *indexed_pixels = arena_alloc(arena, pixel_count * sizeof(GifColorIndex));
    if (indexed_pixels == NULL) {
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    IndexedColor *indexed_colors = arena_alloc(arena, color_count * sizeof(IndexedColor));
    if (indexed_colors == NULL) {
        return NULL;
    }
    {
        f32 const *color_iter = colors;
        GifColorIndex index = 0;
        while (color_iter != colors + color_count * COMPONENTS_PER_COLOR) {
            indexed_colors[index].color_index = index;
            indexed_colors[index].color = (f32x3){color_iter[0], color_iter[1], color_iter[2]};

            index += 1;
            color_iter += COMPONENTS_PER_COLOR;
        }
    }

    if (arena_memory_left(arena) < color_count * sizeof(ColorTreeNode)) {
        return NULL;
    }
    ColorTreeNode *color_tree = color_tree_construct(indexed_colors, color_count, 0, arena);

    f32 const *pixel_iter = pixels;
    GifColorIndex *indexed_pixel_iter = indexed_pixels;
    while (pixel_iter != pixels + pixel_count * COMPONENTS_PER_COLOR) {
        GifColorIndex nearest_color_index;
        f32 distance_to_nearest = INFINITY;
        color_tree_get_nearest(
            color_tree,
            (f32x3){pixel_iter[0], pixel_iter[1], pixel_iter[2]},
            0,
            &distance_to_nearest,
            &nearest_color_index
        );
        *indexed_pixel_iter = nearest_color_index;

        pixel_iter += COMPONENTS_PER_COLOR;
        indexed_pixel_iter += 1;
    }

    arena_rewind(arena, snapshot);
    return indexed_pixels;
}

GifOutputBuffer *gif_out_buffer_create(isize min_capacity, void *arena) {
    isize capacity = min_capacity;
    if (capacity < GIF_OUT_BUFFER_MIN_CAPACITY) {
        capacity = GIF_OUT_BUFFER_MIN_CAPACITY;
    }

    GifOutputBuffer *out_buffer = arena_alloc(arena, sizeof(GifOutputBuffer));
    if (out_buffer == NULL) {
        return NULL;
    }

    out_buffer->data = arena_alloc(arena, capacity);
    if (out_buffer->data == NULL) {
        return NULL;
    }
    memset(out_buffer->data, 0, (size_t)capacity);

    out_buffer->encoded_size = 0;
    out_buffer->byte_pos = 0;
    out_buffer->bit_pos = 0;
    out_buffer->capacity = capacity;

    return out_buffer;
}

isize gif_out_buffer_capacity_left(GifOutputBuffer const *out_buffer) {
    return out_buffer->capacity - out_buffer->byte_pos - (out_buffer->bit_pos == 0 ? 0 : 1);
}

bool gif_out_buffer_grow(GifOutputBuffer *out_buffer, isize min_capacity, void *arena) {
    isize new_capacity = isize_max(
        min_capacity,
        out_buffer->capacity + GIF_OUT_BUFFER_MIN_CAPACITY
    );
    assert(new_capacity > out_buffer->capacity);

    u8 *new_buffer = arena_realloc(
        arena,
        out_buffer->data,
        out_buffer->capacity,
        new_capacity
    );
    if (new_buffer == NULL) {
        return false;
    }
    memset(new_buffer + out_buffer->capacity, 0, (size_t)(new_capacity - out_buffer->capacity));

    out_buffer->data = new_buffer;
    out_buffer->capacity = new_capacity;

    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    return true;
}

void gif_out_buffer_reset(GifOutputBuffer *out_buffer) {
    isize unencoded_size =
        (out_buffer->byte_pos - out_buffer->encoded_size) +
        (out_buffer->bit_pos == 0 ? 0 : 1);
    memmove(out_buffer->data, out_buffer->data + out_buffer->encoded_size, (size_t)unencoded_size);

    out_buffer->byte_pos -= out_buffer->encoded_size;
    out_buffer->encoded_size = 0;

    u8 *free_memory_begin =
        out_buffer->data + out_buffer->byte_pos + (out_buffer->bit_pos == 0 ? 0 : 1);
    u8 *free_memory_end = out_buffer->data + out_buffer->capacity;
    memset(free_memory_begin, 0, (size_t)(free_memory_end - free_memory_begin));
}

typedef struct {
    GifColorIndex *indices;
    isize count;
} ColorSequence;

typedef struct {
    GifColorIndex *data;
    isize size;
    isize capacity;
} ColorSequenceArena;

static GifColorIndex *color_sequence_alloc(ColorSequenceArena *arena, isize size) {
    assert(arena->size + size < arena->capacity);

    GifColorIndex *indices = arena->data + arena->size;
    arena->size += size;
    return indices;
}

// Based on djb2 hashing function:
// http://www.cse.yorku.ca/~oz/hash.html
static u32 color_sequence_hash(ColorSequence sequence) {
    u32 hash = 5381;

    GifColorIndex *sequence_iter = sequence.indices;
    while (sequence_iter < sequence.indices + sequence.count) {
        hash = ((hash << 5) + hash) + *sequence_iter;

        sequence_iter += 1;
    }

    return hash;
}

static bool color_sequence_equals(ColorSequence left, ColorSequence right) {
    assert(left.count != 0 || right.count != 0);

    if (left.count != right.count) {
        return false;
    }

    return memcmp(left.indices, right.indices, (size_t)(left.count * sizeof(GifColorIndex))) == 0;
}

typedef struct {
    GifColorIndex *indices;
    isize count;
    isize capacity;
} ColorArray;

void color_array_push(ColorArray *colors, GifColorIndex index) {
    assert(colors->count < colors->capacity);

    colors->indices[colors->count] = index;
    colors->count += 1;
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
    LzwTableBucket *buckets;
    isize sequence_count;
    isize bucket_count;
} LzwTable;

static void lzw_table_insert(LzwTable *lzw_table, ColorSequence new_sequence, LzwCode new_code) {
    assert(new_sequence.count > 0);
    assert(lzw_table->sequence_count < lzw_table->bucket_count);

    isize bucket_index = color_sequence_hash(new_sequence) % lzw_table->bucket_count;
    while (lzw_table->buckets[bucket_index].sequence.count > 0) {
        assert(!color_sequence_equals(new_sequence, lzw_table->buckets[bucket_index].sequence));

        bucket_index = (bucket_index + 1) % lzw_table->bucket_count;
    }

    lzw_table->buckets[bucket_index].sequence = new_sequence;
    lzw_table->buckets[bucket_index].code = new_code;

    lzw_table->sequence_count += 1;
}

static bool lzw_table_get_code(
    LzwTable const *lzw_table,
    ColorSequence needle_sequence,
    LzwCode *code
) {
    isize bucket_index = color_sequence_hash(needle_sequence) % lzw_table->bucket_count;
    while (
        lzw_table->buckets[bucket_index].sequence.count > 0 &&
        !color_sequence_equals(needle_sequence, lzw_table->buckets[bucket_index].sequence)
    ) {
        bucket_index = (bucket_index + 1) % lzw_table->bucket_count;
    }

    if (lzw_table->buckets[bucket_index].sequence.count > 0) {
        *code = lzw_table->buckets[bucket_index].code;
        return true;
    } else {
        return false;
    }
}

typedef enum {
    GIF_ENCODER_IDLE,
    GIF_ENCODER_READY_FOR_NEXT_FRAME,
    GIF_ENCODER_FRAME_STARTED,
} GifEncoderState;

struct GifEncoder {
    GifEncoderState state;

    isize width;
    isize height;
    u16 frame_delay;

    // These are the actual numbers of colors provided by the user, not the rounded amounts.
    isize global_color_count;
    isize local_color_count;

    int lzw_code_min_bit_length;
    LzwCode lzw_clear_code;
    LzwCode lzw_end_code;
    int starting_lzw_code_bit_length;
    LzwCode starting_lzw_code;

    // Index into the GIF output buffer relative to (data + encoded_size), so that it does not get
    // invalidated, when the output buffer is getting flushed.
    isize data_block_begin;

    int lzw_code_bit_length;
    // The type is isize and not LzwCode, because in theory this variable could overflow the LzwCode
    // type, when it gets incremented when a new code is inserted. I know it does not matter here,
    // because LzwCode is literally an alias for isize, but still.
    isize next_available_lzw_code;

    ColorSequenceArena color_sequence_arena;
    LzwTable lzw_table;
    ColorArray current_sequence;
};

static void gif_encoder_init(GifEncoder *encoder, Arena *arena) {
    encoder->state = GIF_ENCODER_IDLE;
    encoder->global_color_count = 0;
    encoder->local_color_count = 0;

    // Memory for the color sequences
    // Worst case scenario: large single-color image with a palette of size 2.
    //
    // Initially LZW table is filled with 2 sequences of length 1. Then we have a range of [4..4096)
    // LZW codes available to be associated with color sequences. While encoding a single-color
    // image we will use longer and longer sequences of the same color (starting from a sequence of
    // size 2) until we reach the maximum LZW code.
    //
    // So, in total we will need:
    // 1 + 1 + 2 + 3 + ... + 4093 = 1 + 1 + (2 + 4093) * 4092 / 2 = 8378372 bytes ~ 8 MiB
    //         \________________/
    //           4092 sequences

    // Minus 2 for the sequences of length 1.
    // Minus 2 for the clear and end codes.
    isize max_sequence_size = 1 + ((MAX_LZW_CODE + 1) - 4);
    isize total_color_count = 1 + 1 + (2 + max_sequence_size) * (max_sequence_size - 2 + 1) / 2;
    u8 *sequence_memory = arena_alloc(arena, total_color_count * sizeof(GifColorIndex));
    encoder->color_sequence_arena.data = sequence_memory;
    encoder->color_sequence_arena.size = 0;
    encoder->color_sequence_arena.capacity = total_color_count;

    // We will need at most 4096 buckets (24 bytes per bucket), allocate twice this amount:
    // 24 * 4096 * 2 = 196608 bytes = 192 KiB
    isize lzw_bucket_count = 2 * (MAX_LZW_CODE + 1);
    encoder->lzw_table.buckets = arena_alloc(arena, lzw_bucket_count * sizeof(LzwTableBucket));
    encoder->lzw_table.sequence_count = 0;
    encoder->lzw_table.bucket_count = lzw_bucket_count;

    // Additionally 4093 bytes for the max sequence of the currently accumulated image colors.
    encoder->current_sequence.indices = arena_alloc(arena, max_sequence_size);
    encoder->current_sequence.count = 0;
    encoder->current_sequence.capacity = max_sequence_size;
}

isize gif_encoder_required_memory(void) {
    uptr const UPTR_MAX = ~(uptr)0;
    Arena fake_arena = {
        .begin = 0,
        .end = (u8 *)UPTR_MAX,
    };

    GifEncoder encoder = {0};
    arena_alloc(&fake_arena, sizeof(GifEncoder));
    gif_encoder_init(&encoder, &fake_arena);

    return (isize)fake_arena.begin;
}

GifEncoder *gif_encoder_create(void *arena) {
    isize required_memory = gif_encoder_required_memory();
    if (arena_memory_left(arena) < required_memory) {
        return NULL;
    }

    GifEncoder *encoder = arena_alloc(arena, sizeof(GifEncoder));
    gif_encoder_init(encoder, arena);

    // Amounts to a delay of 1 second.
    encoder->frame_delay = 100;

    return encoder;
}

static inline isize u8_write(u8 value, u8 *dest) {
    dest[0] = value;
    return 1;
}

static inline isize u16_write_le(u16 value, u8 *dest) {
    dest[0] = (u8)(value & 0xff);
    dest[1] = (u8)(value >> 8);
    return 2;
}

static inline int u32_leading_zeros(u32 value) {
#if defined(__clang__) || defined(__GNUC__)
    return __builtin_clz(value);
#elif defined(_MSC_VER)
    return (int)__lzcnt(value);
#else
    int result = 0;
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

static inline int u32_log2_ceil(u32 value) {
    assert(value != 0);
    return 32 - u32_leading_zeros(value - 1);
}

static inline bool isize_power_of_two(isize value) {
    assert(value > 0);
    return (value & (value - 1)) == 0;
}

// GIF color tables can only have sizes which are powers of two >= 2.
static inline isize gif_color_count_round_up(isize color_count) {
    if (color_count == 0) {
        return 0;
    } else if (color_count <= 2) {
        // From Appendix F of GIF89a specification:
        // The first byte of the Compressed Data stream is a value indicating the minimum number of
        // bits required to represent the set of actual pixel values. Normally this will be the same
        // as the number of color bits. Because of some algorithmic constraints however, black &
        // white images which have one color bit must be indicated as having a code size of 2.
        //
        // 4 because log₂(4) is equal to 2.
        return 4;
    } else {
        return 1 << u32_log2_ceil((u32)color_count);
    }
}

void gif_encoder_start(
    GifEncoder *encoder,
    isize width,
    isize height,
    u8 const *global_colors,
    isize global_color_count,
    GifOutputBuffer *out_buffer
) {
    assert(encoder->state == GIF_ENCODER_IDLE);
    assert(width <= GIF_MAX_WIDTH && height <= GIF_MAX_HEIGHT);
    assert(global_color_count == 0 || global_colors != NULL);
    assert(global_color_count <= GIF_MAX_COLORS);

    isize global_color_count_rounded = gif_color_count_round_up(global_color_count);

    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);
    assert(out_buffer->bit_pos == 0);
    u8 *out_iter = out_buffer->data + out_buffer->byte_pos;

    // GIF header
    // -----------------
    // 6 bytes          for the "GIF89a" header
    // 7 bytes          for a logical screen descriptor
    // 256 * 3 bytes    for a max size global color table
    // 19 bytes         for the Netscape Looping Application Extension
    // -----------------
    // 800 bytes

    // GIF89a magic

    out_iter += u8_write('G', out_iter);
    out_iter += u8_write('I', out_iter);
    out_iter += u8_write('F', out_iter);
    out_iter += u8_write('8', out_iter);
    out_iter += u8_write('9', out_iter);
    out_iter += u8_write('a', out_iter);

    // Logical screen descriptor

    out_iter += u16_write_le((u16)width, out_iter);
    out_iter += u16_write_le((u16)height, out_iter);

    // |A|BBB|C|DDD|
    // A: Whether global color table is present
    // D: Size of the global color table is: 2^(D+1)
    if (global_color_count > 0) {
        out_iter += u8_write(
            0x80 | (u8)(u32_log2_ceil((u32)global_color_count_rounded) - 1),
            out_iter
        );
    } else {
        out_iter += u8_write(0x00, out_iter);
    }

    // Background color index (from the global color table)
    out_iter += u8_write(0, out_iter);
    // Unimportant
    out_iter += u8_write(0, out_iter);

    // Global color table

    u8 const *global_color_iter = global_colors;
    for (isize i = 0; i < global_color_count_rounded; i += 1) {
        if (global_color_iter < global_colors + global_color_count * COMPONENTS_PER_COLOR) {
            out_iter += u8_write(global_color_iter[0], out_iter);
            out_iter += u8_write(global_color_iter[1], out_iter);
            out_iter += u8_write(global_color_iter[2], out_iter);

            global_color_iter += COMPONENTS_PER_COLOR;
        } else {
            out_iter += u8_write(0x00, out_iter);
            out_iter += u8_write(0x00, out_iter);
            out_iter += u8_write(0x00, out_iter);
        }
    }

    // Netscape Looping Application Extension
    // https://web.archive.org/web/20180620143135/https://www.vurdalakov.net/misc/gif/netscape-looping-application-extension

    // Extension introducer
    out_iter += u8_write(0x21, out_iter);
    // Application extension label
    out_iter += u8_write(0xff, out_iter);
    // Extension size + data
    {
        char extension_data[] = "NETSCAPE2.0";
        out_iter += u8_write(lengthof(extension_data), out_iter);
        memcpy(out_iter, extension_data, 11);
        out_iter += lengthof(extension_data);
    }
    // Data sub-blocks
    out_iter += u8_write(3, out_iter);
    out_iter += u8_write(0x01, out_iter);
    // Loop count (0 for infinity)
    out_iter += u16_write_le(0, out_iter);
    out_iter += u8_write(0, out_iter);

    // Update the output buffer

    assert(out_iter <= out_buffer->data + out_buffer->capacity);
    out_buffer->encoded_size = out_iter - out_buffer->data;
    out_buffer->byte_pos = out_iter - out_buffer->data;

    // Update encoder state

    encoder->state = GIF_ENCODER_READY_FOR_NEXT_FRAME;
    encoder->width = width;
    encoder->height = height;
    encoder->global_color_count = global_color_count;
}

void gif_frame_delay(GifEncoder *encoder, f32 seconds) {
    encoder->frame_delay = (u16)(roundf(seconds * 100.0F));
    if (encoder->frame_delay == 0) {
        encoder->frame_delay = 1;
    }
}

static void gif_encoder_write_lzw_code(
    GifEncoder *encoder,
    LzwCode code,
    int code_bit_length,
    GifOutputBuffer *out_buffer
) {
    assert(encoder->state == GIF_ENCODER_FRAME_STARTED);
    assert(code_bit_length <= LZW_CODE_MAX_BIT_LENGTH);
    assert(code <= MAX_LZW_CODE);
    {
        isize bits_left = (out_buffer->capacity - out_buffer->byte_pos) * 8 - out_buffer->bit_pos;
        assert(bits_left >= code_bit_length);
    }

    if (out_buffer->byte_pos - (out_buffer->encoded_size + encoder->data_block_begin) == 256) {
        assert(out_buffer->bit_pos == 0);
        assert(out_buffer->capacity - out_buffer->byte_pos >= 1 + (code_bit_length + 7) / 8);

        out_buffer->data[out_buffer->encoded_size + encoder->data_block_begin] = 255;
        out_buffer->encoded_size += 256;
        out_buffer->byte_pos += 1;
        encoder->data_block_begin = 0;
    }

    int code_bits_left = code_bit_length;
    LzwCode code_remainder = code;

    out_buffer->data[out_buffer->byte_pos] |= (u8)((code_remainder << out_buffer->bit_pos) & 0xff);

    if (code_bits_left + out_buffer->bit_pos < 8) {
        out_buffer->bit_pos += code_bits_left;
        code_bits_left = 0;
    } else {
        isize bits_written = 8 - out_buffer->bit_pos;

        out_buffer->bit_pos = 0;
        code_bits_left -= bits_written;
        code_remainder >>= bits_written;

        out_buffer->byte_pos += 1;
    }

    while (code_bits_left > 0) {
        assert(out_buffer->bit_pos == 0);

        if (out_buffer->byte_pos - (out_buffer->encoded_size + encoder->data_block_begin) == 256) {
            assert(out_buffer->capacity - out_buffer->byte_pos >= 1 + (code_bits_left + 7) / 8);

            out_buffer->data[out_buffer->encoded_size + encoder->data_block_begin] = 255;
            out_buffer->encoded_size += 256;
            out_buffer->byte_pos += 1;
            encoder->data_block_begin = 0;
        }

        out_buffer->data[out_buffer->byte_pos] |= (u8)(code_remainder & 0xff);

        if (code_bits_left < 8) {
            out_buffer->bit_pos = code_bits_left;
            code_bits_left = 0;
        } else {
            isize bits_written = 8;

            out_buffer->bit_pos = 0;
            code_bits_left -= bits_written;
            code_remainder >>= bits_written;

            out_buffer->byte_pos += 1;
        }
    }
}

static void gif_encoder_reset_lzw_table(GifEncoder *encoder) {
    encoder->color_sequence_arena.size = 0;
    encoder->lzw_table.sequence_count = 0;
    memset(
        encoder->lzw_table.buckets,
        0,
        (size_t)(encoder->lzw_table.bucket_count * sizeof(LzwTableBucket))
    );
    for (
        isize color_index = 0;
        color_index < (1 << encoder->lzw_code_min_bit_length);
        color_index += 1
    ) {
        ColorSequence single_color_sequence = {
            .count = 1,
            .indices = color_sequence_alloc(&encoder->color_sequence_arena, 1),
        };
        single_color_sequence.indices[0] = (GifColorIndex)color_index;

        lzw_table_insert(&encoder->lzw_table, single_color_sequence, color_index);
    }
}

void gif_encoder_start_frame(
    GifEncoder *encoder,
    u8 const *local_colors,
    isize local_color_count,
    GifOutputBuffer *out_buffer
) {
    assert(encoder->state == GIF_ENCODER_READY_FOR_NEXT_FRAME);
    assert(local_color_count == 0 || local_colors != NULL);
    assert(local_color_count <= GIF_MAX_COLORS);
    assert(encoder->global_color_count > 0 || local_color_count > 0);

    isize local_color_count_rounded = gif_color_count_round_up(local_color_count);

    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);
    assert(out_buffer->bit_pos == 0);
    u8 *out_iter = out_buffer->data + out_buffer->byte_pos;

    // Image header
    // -----------------
    // 8 bytes          for a graphic control extension
    // 10 bytes         for an image descriptor
    // 256 * 3 bytes    for a max size local color table
    // 1 byte           for a minimum LZW code size
    // 1 byte           for the size of the first image data block
    // 9 bits           for the LZW clear code (9 bits in case of max size color table)
    // -----------------
    // 790 bytes

    // Graphic control extension

    // Extension introducer
    out_iter += u8_write(0x21, out_iter);
    // Graphic control extension label
    out_iter += u8_write(0xf9, out_iter);
    // Extension data size excluding the terminator
    out_iter += u8_write(4, out_iter);
    // |AAA|BBB|C|D|
    // B: "2" means clear the canvas with background color after each frame
    // D: Whether transparency is enabled
    out_iter += u8_write(0x08, out_iter);
    // Frame delay measured in 1/100ths of a second
    out_iter += u16_write_le(encoder->frame_delay, out_iter);
    // Transparent color index (only matters if transparency is enabled)
    out_iter += u8_write(0, out_iter);
    // Terminator (a zero-size data block)
    out_iter += u8_write(0, out_iter);

    // Image descriptor

    // Image descriptor introducer
    out_iter += u8_write(0x2c, out_iter);
    // Image left position
    out_iter += u16_write_le(0, out_iter);
    // Image top position
    out_iter += u16_write_le(0, out_iter);

    out_iter += u16_write_le((u16)encoder->width, out_iter);
    out_iter += u16_write_le((u16)encoder->height, out_iter);

    // |A|B|C|DD|EEE|
    // A: Whether a local color table is present
    // E: Size of the local color table is: 2^(E+1)
    if (local_color_count > 0) {
        out_iter += u8_write(
            0x80 | (u8)(u32_log2_ceil((u32)local_color_count_rounded) - 1),
            out_iter
        );
    } else {
        out_iter += u8_write(0x00, out_iter);
    }

    // Local color table

    u8 const *local_color_iter = local_colors;
    for (isize i = 0; i < local_color_count_rounded; i += 1) {
        if (local_color_iter < local_colors + local_color_count * COMPONENTS_PER_COLOR) {
            out_iter += u8_write(local_color_iter[0], out_iter);
            out_iter += u8_write(local_color_iter[1], out_iter);
            out_iter += u8_write(local_color_iter[2], out_iter);

            local_color_iter += COMPONENTS_PER_COLOR;
        } else {
            out_iter += u8_write(0x00, out_iter);
            out_iter += u8_write(0x00, out_iter);
            out_iter += u8_write(0x00, out_iter);
        }
    }

    // Minimum LZW code bit length

    int lzw_code_min_bit_length;
    if (local_color_count > 0) {
        lzw_code_min_bit_length = u32_log2_ceil((u32)local_color_count_rounded);
    } else {
        isize global_color_count_rounded = gif_color_count_round_up(encoder->global_color_count);
        lzw_code_min_bit_length = u32_log2_ceil((u32)global_color_count_rounded);
    }
    out_iter += u8_write((u8)lzw_code_min_bit_length, out_iter);

    // Size of the first image data block
    // Just skip over it for now, because we don't know the size of the data block upfront.
    out_iter += 1;

    // Update encoder state

    encoder->state = GIF_ENCODER_FRAME_STARTED;
    encoder->local_color_count = local_color_count;

    encoder->lzw_code_min_bit_length = lzw_code_min_bit_length;
    encoder->lzw_clear_code = 1 << lzw_code_min_bit_length;
    encoder->lzw_end_code = encoder->lzw_clear_code + 1;
    encoder->starting_lzw_code = encoder->lzw_clear_code + 2;
    encoder->starting_lzw_code_bit_length = lzw_code_min_bit_length + 1;

    encoder->lzw_code_bit_length = encoder->starting_lzw_code_bit_length;
    encoder->next_available_lzw_code = encoder->starting_lzw_code;

    assert(out_iter <= out_buffer->data + out_buffer->capacity);
    // Subtract one to exclude the data block size (which is now unknown) from the encoded data.
    out_buffer->encoded_size = out_iter - out_buffer->data - 1;
    out_buffer->byte_pos = out_iter - out_buffer->data;
    encoder->data_block_begin = 0;

    encoder->current_sequence.count = 0;

    gif_encoder_reset_lzw_table(encoder);

    // Output a clear code because the spec says so:
    // > encoders should output a clear code as the first code of each image data stream.

    gif_encoder_write_lzw_code(
        encoder,
        encoder->lzw_clear_code,
        encoder->starting_lzw_code_bit_length,
        out_buffer
    );
}

isize gif_encoder_feed_frame(
    GifEncoder *encoder,
    GifColorIndex const *pixels,
    isize pixel_count,
    GifOutputBuffer *out_buffer
) {
    assert(encoder->state == GIF_ENCODER_FRAME_STARTED);
    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    // The user must indicate the end of the frame explicitly.
    if (pixel_count == 0) {
        return 0;
    }

    assert(encoder->global_color_count > 0 || encoder->local_color_count > 0);
    isize color_count = encoder->global_color_count;
    if (encoder->local_color_count > 0) {
        color_count = encoder->local_color_count;
    }

    GifColorIndex const *pixel_iter = pixels;
    GifColorIndex const *pixels_end = pixels + pixel_count;

    // This should only happen once, right after the frame has been started.
    //
    // Cannot do this inside of gif_encoder_start_frame, because at that point the image pixels are
    // not available yet.
    if (encoder->current_sequence.count == 0) {
        assert(*pixel_iter < color_count);
        color_array_push(&encoder->current_sequence, *pixel_iter);
        pixel_iter += 1;
    }

    while (pixel_iter < pixels_end) {
        // Make sure that we will be able to exit this function without a buffer overflow.
        //
        // Conservative check:
        //
        // + (encoder->lzw_code_bit_length) bits
        // just to write the next code itself.
        //
        // + 1 byte for the data block size
        // in case we overflow the current data block and we need to start a next one.
        //
        // + LZW_CODE_MAX_BIT_LENGTH bits for the clear code
        // in case we reach the maximum LZW code

        isize bits_left = (out_buffer->capacity - out_buffer->byte_pos) * 8 - out_buffer->bit_pos;
        if (bits_left <= encoder->lzw_code_bit_length + 8 + LZW_CODE_MAX_BIT_LENGTH) {
            break;
        }

        // It should get initialized on the first iteration of the loop, because the LZW table
        // contains every color sequence of length 1.
        LzwCode next_lzw_code;
        bool new_sequence_found = !lzw_table_get_code(
            &encoder->lzw_table,
            (ColorSequence){encoder->current_sequence.indices, encoder->current_sequence.count},
            &next_lzw_code
        );
        // 1st case: This is the first call to the gif_encoder_feed_frame and we explicitly
        // initialized current_sequence to the first image pixel. Sequences of size 1 are guaranteed
        // to exist in the LZW table.
        //
        // 2nd case: These are pixels left unencoded after the last call to the
        // gif_encoder_feed_frame, which means that they are guaranteed to exist in the LZW table
        // (otherwise we would have had put them into the table).
        assert(!new_sequence_found);

        // Accumulate pixels from the image until we find a sequence which is not yet in the table.
        do {
            assert(*pixel_iter < color_count);
            color_array_push(&encoder->current_sequence, *pixel_iter);
            pixel_iter += 1;

            new_sequence_found = !lzw_table_get_code(
                &encoder->lzw_table,
                (ColorSequence){encoder->current_sequence.indices, encoder->current_sequence.count},
                &next_lzw_code
            );
        } while (pixel_iter < pixels_end && !new_sequence_found);

        if (new_sequence_found) {
            gif_encoder_write_lzw_code(
                encoder,
                next_lzw_code,
                encoder->lzw_code_bit_length,
                out_buffer
            );

            // All possible sequences of length 1 are guaranteed to already exist in the table.
            assert(encoder->current_sequence.count >= 2);
            assert(encoder->next_available_lzw_code <= MAX_LZW_CODE);

            GifColorIndex *new_sequence_indices = color_sequence_alloc(
                &encoder->color_sequence_arena,
                encoder->current_sequence.count
            );
            ColorSequence new_sequence = {
                .count = encoder->current_sequence.count,
                .indices = new_sequence_indices,
            };
            memmove(
                new_sequence.indices,
                encoder->current_sequence.indices,
                (size_t)(encoder->current_sequence.count * sizeof(GifColorIndex))
            );

            // Increase the LZW code length here because the spec says so:
            // > Whenever the LZW code value would exceed the current code length,
            // > the code length is increased by one. The packing/unpacking of these
            // > codes must then be altered to reflect the new code length.

            if (isize_power_of_two(encoder->next_available_lzw_code)) {
                encoder->lzw_code_bit_length += 1;
            }

            lzw_table_insert(
                &encoder->lzw_table,
                new_sequence,
                encoder->next_available_lzw_code
            );
            encoder->next_available_lzw_code += 1;

            GifColorIndex next_after_existing =
                encoder->current_sequence.indices[encoder->current_sequence.count - 1];
            encoder->current_sequence.count = 1;
            encoder->current_sequence.indices[0] = next_after_existing;

            // When the table is full, the encoder can chose to use the table as is, making no
            // changes to it until the encoder chooses to clear it. The encoder during this time
            // sends out codes that are of the maximum Code Size.
            //
            // Which means that you don't have to necessarily emit the clear code here, but clearing
            // the LZW table usually leads to better compression.
            if (encoder->next_available_lzw_code > MAX_LZW_CODE) {
                gif_encoder_write_lzw_code(
                    encoder,
                    encoder->lzw_clear_code,
                    encoder->lzw_code_bit_length,
                    out_buffer
                );

                // Dealloc color sequences and empty the LZW table
                gif_encoder_reset_lzw_table(encoder);
                encoder->lzw_code_bit_length = encoder->starting_lzw_code_bit_length;
                encoder->next_available_lzw_code = encoder->starting_lzw_code;
            }
        }
    }

    return pixel_iter - pixels;
}

void gif_encoder_finish_frame(GifEncoder *encoder, GifOutputBuffer *out_buffer) {
    assert(encoder->state == GIF_ENCODER_FRAME_STARTED);
    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    // Encode the leftover sequence.
    if (encoder->current_sequence.count > 0) {
        LzwCode last_lzw_code;
        bool last_sequence_exists = lzw_table_get_code(
            &encoder->lzw_table,
            (ColorSequence){encoder->current_sequence.indices, encoder->current_sequence.count},
            &last_lzw_code
        );
        assert(last_sequence_exists);
        gif_encoder_write_lzw_code(
            encoder,
            last_lzw_code,
            encoder->lzw_code_bit_length,
            out_buffer
        );
        encoder->current_sequence.count = 0;
    }

    gif_encoder_write_lzw_code(
        encoder,
        encoder->lzw_end_code,
        encoder->lzw_code_bit_length,
        out_buffer
    );
    if (out_buffer->bit_pos != 0) {
        out_buffer->bit_pos = 0;
        out_buffer->byte_pos += 1;
    }

    isize data_block_pos = out_buffer->encoded_size + encoder->data_block_begin;
    // -1 to exclude the size byte from the total size of the block.
    out_buffer->data[data_block_pos] = (u8)(out_buffer->byte_pos - data_block_pos - 1);

    // Terminator block
    out_buffer->data[out_buffer->byte_pos] = 0x00;
    out_buffer->byte_pos += 1;

    out_buffer->encoded_size = out_buffer->byte_pos;

    encoder->state = GIF_ENCODER_READY_FOR_NEXT_FRAME;
}

bool gif_encode_whole_frame(
    GifEncoder *encoder,
    u8 const *local_colors,
    isize local_color_count,
    GifColorIndex *pixels,
    GifOutputBuffer *out_buffer
) {
    assert(encoder->state == GIF_ENCODER_READY_FOR_NEXT_FRAME);

    GifOutputBuffer out_buffer_rewind = *out_buffer;

    GifColorIndex *pixel_iter = pixels;
    GifColorIndex *pixels_end = pixels + encoder->width * encoder->height;

    if (gif_out_buffer_capacity_left(out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
        *out_buffer = out_buffer_rewind;
        return false;
    }
    gif_encoder_start_frame(encoder, local_colors, local_color_count, out_buffer);

    while (pixel_iter != pixels_end) {
        if (gif_out_buffer_capacity_left(out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
            *out_buffer = out_buffer_rewind;
            return false;
        }
        pixel_iter += gif_encoder_feed_frame(
            encoder,
            pixel_iter,
            pixels_end - pixel_iter,
            out_buffer
        );
    }

    if (gif_out_buffer_capacity_left(out_buffer) < GIF_OUT_BUFFER_MIN_CAPACITY) {
        *out_buffer = out_buffer_rewind;
        return false;
    }
    gif_encoder_finish_frame(encoder, out_buffer);

    return true;
}

void gif_encoder_finish(GifEncoder *encoder, GifOutputBuffer *out_buffer) {
    assert(encoder->state == GIF_ENCODER_READY_FOR_NEXT_FRAME);
    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    // GIF trailer
    out_buffer->data[out_buffer->byte_pos] = 0x3b;
    out_buffer->byte_pos += 1;

    out_buffer->encoded_size = out_buffer->byte_pos;

    encoder->state = GIF_ENCODER_IDLE;
}
