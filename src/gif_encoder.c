// GIF references:
// - https://www.w3.org/Graphics/GIF/spec-gif89a.txt
// - https://en.wikipedia.org/wiki/GIF
// - https://web.archive.org/web/20180620143135/https://www.vurdalakov.net/misc/gif/netscape-looping-application-extension

// TODO: Transdiff from ffmpeg: unchanged pixels are encoded as transparent on next frame.
// TODO: Color search: sort colors along 3 axis, first check colors which have closer projections.
// TODO: Pass the number of color components in arguments (1 for monochrome, 3 for RGB, 4 for RGBA).

#include "gif_encoder.h"

#include <assert.h> // assert
#include <stdlib.h> // qsort
#include <stddef.h> // NULL
#include <math.h>   // copysignf, powf, roundf, cbrtf, INFINITY
#include <string.h> // memmove, memcpy, memset

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


#define sizeof(expr) (isize)sizeof(expr)
#define lengthof(string) (sizeof(string) - 1)
#define countof(expr) (sizeof(expr) / sizeof((expr)[0]))

#define USIZE_BITS ((int)sizeof(usize) * 8)

static inline isize u8_write(u8 value, u8 *dest) {
    dest[0] = value;
    return 1;
}

static inline isize u16_write_le(u16 value, u8 *dest) {
    dest[0] = (u8)(value & 0xff);
    dest[1] = (u8)(value >> 8);
    return 2;
}

#if defined(__clang__) || defined(__GNUC__)
    #define u32_leading_zeroes __builtin_clz
    #define u64_leading_zeroes __builtin_clzll
#elif defined(_MSC_VER)
    #define u32_leading_zeroes __lzcnt
    #define u64_leading_zeroes __lzcnt64
#else
    #define u32_leading_zeroes u32_leading_zeroes_impl
    #define u64_leading_zeroes u64_leading_zeroes_impl
#endif

static int u32_leading_zeroes_impl(u32 value) {
    int leading_zeroes = 0;
    while (value != 0) {
        value >>= 1;
        leading_zeroes += 1;
    }
    return 32 - leading_zeroes;
}

static int u64_leading_zeroes_impl(u64 value) {
    int leading_zeroes = 0;
    while (value != 0) {
        value >>= 1;
        leading_zeroes += 1;
    }
    return 64 - leading_zeroes;
}

static inline int u32_log2_ceil(u32 value) {
    assert(value != 0);
    return 32 - u32_leading_zeroes(value - 1);
}

static inline int usize_leading_zeroes(usize value) {
    if (sizeof(usize) == 8) {
        return u64_leading_zeroes(value);
    } else {
        return u32_leading_zeroes(value);
    }
}

static inline isize isize_min(isize left, isize right) {
    return left < right ? left : right;
}

static inline isize isize_max(isize left, isize right) {
    return left > right ? left : right;
}

static inline isize isize_clamp(isize value, isize min, isize max) {
    if (value < min) {
        return min;
    }
    if (value > max) {
        return max;
    }
    return value;
}

static inline isize isize_abs(isize value) {
    if (value >= 0) {
        return value;
    } else {
        return -value;
    }
}

static inline bool isize_power_of_two(isize value) {
    assert(value > 0);
    return (value & (value - 1)) == 0;
}

static inline f32 f32_min(f32 left, f32 right) {
    return left < right ? left : right;
}

static inline f32 f32_max(f32 left, f32 right) {
    return left > right ? left : right;
}

typedef struct {
    f32 x;
    f32 y;
    f32 z;
} f32x3;

typedef struct {
    f32x3 min;
    f32x3 max;
} f32xbox;

static inline f32x3 f32x3_max(f32x3 left, f32x3 right) {
    return (f32x3){
        left.x > right.x ? left.x : right.x,
        left.y > right.y ? left.y : right.y,
        left.z > right.z ? left.z : right.z,
    };
}

static inline f32x3 f32x3_min(f32x3 left, f32x3 right) {
    return (f32x3){
        left.x < right.x ? left.x : right.x,
        left.y < right.y ? left.y : right.y,
        left.z < right.z ? left.z : right.z,
    };
}

static inline f32x3 f32x3_add(f32x3 left, f32x3 right) {
    return (f32x3){left.x + right.x, left.y + right.y, left.z + right.z};
}

static inline f32x3 f32x3_sub(f32x3 left, f32x3 right) {
    return (f32x3){left.x - right.x, left.y - right.y, left.z - right.z};
}

static inline f32x3 f32x3_mul(f32x3 left, f32x3 right) {
    return (f32x3){left.x * right.x, left.y * right.y, left.z * right.z};
}

static inline f32 f32x3_dot(f32x3 left, f32x3 right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

static inline f32x3 f32x3_scale(f32x3 vector, f32 scalar) {
    return (f32x3){vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

static inline f32x3 f32x3_clamp(f32x3 value, f32x3 min, f32x3 max) {
    if (value.x < min.x) { value.x = min.x; }
    if (value.y < min.y) { value.y = min.y; }
    if (value.z < min.z) { value.z = min.z; }

    if (value.x > max.x) { value.x = max.x; }
    if (value.y > max.y) { value.y = max.y; }
    if (value.z > max.z) { value.z = max.z; }

    return value;
}

static int f32x3_compare_by_x(void const *left, void const *right) {
    return (int)copysignf(1.0F, ((f32x3 *)left)->x - ((f32x3 *)right)->x);
}

static int f32x3_compare_by_y(void const *left, void const *right) {
    return (int)copysignf(1.0F, ((f32x3 *)left)->y - ((f32x3 *)right)->y);
}

static int f32x3_compare_by_z(void const *left, void const *right) {
    return (int)copysignf(1.0F, ((f32x3 *)left)->z - ((f32x3 *)right)->z);
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
    assert(((uptr)old_ptr & (ARENA_ALIGNMENT - 1)) == 0);

    if (new_size == 0) {
        return NULL;
    }
    if (old_size >= new_size) {
        return old_ptr;
    }
    if (old_ptr == NULL) {
        return arena_alloc(arena, new_size);
    }

    if ((u8 *)old_ptr + old_size == arena->begin) {
        arena->begin = old_ptr;
        return arena_alloc(arena, new_size);
    } else {
        void *new_ptr = arena_alloc(arena, new_size);
        if (new_ptr == NULL) {
            return NULL;
        }

        memcpy(new_ptr, old_ptr, old_size);
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

    u8 const *srgb_color_iter = srgb_colors;
    f32 *float_color_iter = float_colors;

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

    f32 const *float_color_iter = float_colors;
    u8 *srgb_color_iter = srgb_colors;

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

    u8 const *srgb_color_iter = srgb_colors;
    f32 *linear_color_iter = linear_colors;

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

    f32 const *linear_color_iter = linear_colors;
    u8 *srgb_color_iter = srgb_colors;

    while (srgb_color_iter != srgb_colors + color_count * COMPONENTS_PER_COLOR) {
        srgb_color_iter[0] = (u8)(linear_component_to_srgb(linear_color_iter[0]) * 255.0F);
        srgb_color_iter[1] = (u8)(linear_component_to_srgb(linear_color_iter[1]) * 255.0F);
        srgb_color_iter[2] = (u8)(linear_component_to_srgb(linear_color_iter[2]) * 255.0F);

        linear_color_iter += COMPONENTS_PER_COLOR;
        srgb_color_iter += COMPONENTS_PER_COLOR;
    }

    return srgb_colors;
}

f32 *srgb_to_lab(u8 const *srgb_colors, isize color_count, void *arena) {
    f32 *lab_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (lab_colors == NULL) {
        return NULL;
    }

    u8 const *srgb_color_iter = srgb_colors;
    f32 *lab_color_iter = lab_colors;

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

    f32 const *lab_color_iter = lab_colors;
    u8 *srgb_color_iter = srgb_colors;

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

    u8 const *srgb_iter = srgb_colors;
    f32 *oklab_iter = oklab_colors;

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

u8 *oklab_to_srgb(f32 const *oklab_colors, isize color_count, void *arena) {
    u8 *srgb_colors = arena_alloc(arena, color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (srgb_colors == NULL) {
        return NULL;
    }

    f32 const *oklab_iter = oklab_colors;
    u8 *srgb_iter = srgb_colors;

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

// Arrays in here are split into chunks. The chunk sizes are as follows: 2^0, 2^1, 2^2, ...
//
// "indices" stores indices into "colors" array and represents a "color hashes -> colors" mapping.
// Since there is going to be at most 2^32 colors, 32-bit indices will do.
// The number of indices always coincides with the capacity, no need to store them separately.
typedef struct {
    i32 *indices[USIZE_BITS];
    isize index_count;

    u8 *colors[USIZE_BITS];
    isize color_count;
    isize color_capacity;
} ColorSet;

// Calculate the location of an item at index "item_index" within the chunked array.
// (e.g. the "item_index"-th index is located at indices[chunk][chunk_pos].)
static inline void chunk_pos(isize item_index, isize *chunk, isize *chunk_pos) {
    *chunk = USIZE_BITS - usize_leading_zeroes(item_index + 1) - 1;
    *chunk_pos = (item_index + 1) - (1 << *chunk);
}

// https://nullprogram.com/blog/2018/07/31/
static u32 color_hash(u8 const *color) {
    u32 hash = (u32)color[0] << 16 | (u32)color[1] << 8 | color[2];

    hash = (hash ^ hash >> 17) * 0xed5ad4bb;
    hash = (hash ^ hash >> 11) * 0xac4c1b51;
    hash = (hash ^ hash >> 15) * 0x31848bab;
    return  hash ^ hash >> 14;
}

// Returns false when it fails to allocate memory from the arena.
static bool color_set_append(
    ColorSet *color_set,
    u8 const *colors, isize color_count,
    void *arena
) {
    u8 const *color_iter = colors;
    u8 const *colors_end = colors + color_count * COMPONENTS_PER_COLOR;
    while (color_iter < colors_end) {
        u8 next_color[3] = {color_iter[0], color_iter[1], color_iter[2]};

        // The load factor has become greater than 2/3.
        if (color_set->color_count * 3 > color_set->index_count * 2) {
            isize new_chunk, new_chunk_pos;
            chunk_pos(color_set->index_count, &new_chunk, &new_chunk_pos);

            isize new_chunk_size = 1 << new_chunk;
            i32 *new_indices = arena_alloc(arena, new_chunk_size * sizeof(i32));
            if (new_indices == NULL) {
                return false;
            }

            color_set->indices[new_chunk] = new_indices;
            color_set->index_count += new_chunk_size;

            for (isize i = 0; i <= new_chunk; i += 1) {
                memset(color_set->indices[i], -1, (1 << i) * sizeof(i32));
            }
            for (i32 color_index = 0; color_index < color_set->color_count; color_index += 1) {
                isize color_chunk, color_chunk_pos;
                chunk_pos(color_index, &color_chunk, &color_chunk_pos);
                u8 *color = &color_set->colors[color_chunk][color_chunk_pos * COMPONENTS_PER_COLOR];

                isize bucket_index = color_hash(color) % color_set->index_count;
                while (true) {
                    isize index_chunk, index_chunk_pos;
                    chunk_pos(bucket_index, &index_chunk, &index_chunk_pos);
                    i32 *bucket = &color_set->indices[index_chunk][index_chunk_pos];
                    if (*bucket < 0) {
                        *bucket = color_index;
                        break;
                    }

                    bucket_index = (bucket_index + 1) % color_set->index_count;
                }
            }
        }
        assert(color_set->color_count * 3 <= color_set->index_count * 2);

        isize bucket_index = color_hash(next_color) % color_set->index_count;
        i32 *bucket;
        while (true) {
            isize index_chunk, index_chunk_pos;
            chunk_pos(bucket_index, &index_chunk, &index_chunk_pos);
            bucket = &color_set->indices[index_chunk][index_chunk_pos];
            if (*bucket < 0) {
                break;
            }

            isize color_chunk, color_chunk_pos;
            chunk_pos(*bucket, &color_chunk, &color_chunk_pos);
            u8 *color = &color_set->colors[color_chunk][color_chunk_pos * COMPONENTS_PER_COLOR];
            if (
                color[0] == next_color[0] && color[1] == next_color[1] && color[2] == next_color[2]
            ) {
                break;
            }

            bucket_index = (bucket_index + 1) % color_set->index_count;
        }

        if (*bucket < 0) {
            // Ran out of the colors array capacity.
            // Only allocate more memory for the colors when we actually need it.
            if (color_set->color_count == color_set->color_capacity) {
                isize new_chunk, new_chunk_pos;
                chunk_pos(color_set->color_capacity, &new_chunk, &new_chunk_pos);
                isize new_chunk_size = 1 << new_chunk;

                u8 *new_colors = arena_alloc(
                    arena,
                    new_chunk_size * COMPONENTS_PER_COLOR * sizeof(u8)
                );
                if (new_colors == NULL) {
                    return false;
                }

                color_set->colors[new_chunk] = new_colors;
                color_set->color_capacity += new_chunk_size;
            }

            isize color_chunk, color_chunk_pos;
            chunk_pos(color_set->color_count, &color_chunk, &color_chunk_pos);
            u8 *color = &color_set->colors[color_chunk][color_chunk_pos * COMPONENTS_PER_COLOR];
            color[0] = next_color[0];
            color[1] = next_color[1];
            color[2] = next_color[2];

            *bucket = (i32)color_set->color_count;

            color_set->color_count += 1;
        }

        color_iter += COMPONENTS_PER_COLOR;
    }

    return true;
}

u8 *colors_unique(u8 const *colors, isize color_count, isize *unique_color_count, void *arena) {
    if (color_count == 0) {
        *unique_color_count = 0;
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    // Must be a positive integer in the for of (2^N - 1).
    isize const INITIAL_CAPACITY = 63;

    // Allocate some memory for the indices and colors from the get-go.
    i32 *initial_indices = arena_alloc(arena, INITIAL_CAPACITY * sizeof(i32));
    u8 *initial_colors = arena_alloc(arena, INITIAL_CAPACITY * COMPONENTS_PER_COLOR * sizeof(u8));
    if (initial_indices == NULL || initial_colors == NULL) {
        *unique_color_count = 0;
        return NULL;
    }
    memset(initial_indices, -1, INITIAL_CAPACITY * sizeof(i32));

    ColorSet color_set = {
        .color_count = 0,
        .color_capacity = INITIAL_CAPACITY,
        .index_count = INITIAL_CAPACITY,
    };
    for (int i = 0; i < USIZE_BITS - usize_leading_zeroes(INITIAL_CAPACITY); i += 1) {
        color_set.colors[i] = &initial_colors[((1 << i) - 1) * COMPONENTS_PER_COLOR];
    }
    for (int i = 0; i < USIZE_BITS - usize_leading_zeroes(INITIAL_CAPACITY); i += 1) {
        color_set.indices[i] = &initial_indices[(1 << i) - 1];
    }

    if (!color_set_append(&color_set, colors, color_count, arena)) {
        *unique_color_count = 0;
        return NULL;
    }

    // Write over the hash set memory.
    arena_rewind(arena, snapshot);
    u8 *unique_colors = arena_alloc(
        arena,
        color_set.color_count * COMPONENTS_PER_COLOR * sizeof(u8)
    );

    u8 *unique_color_iter = unique_colors;
    for (isize color_index = 0; color_index < color_set.color_count; color_index += 1) {
        isize color_chunk, color_chunk_pos;
        chunk_pos(color_index, &color_chunk, &color_chunk_pos);
        u8 *color = &color_set.colors[color_chunk][color_chunk_pos * COMPONENTS_PER_COLOR];

        unique_color_iter[0] = color[0];
        unique_color_iter[1] = color[1];
        unique_color_iter[2] = color[2];

        unique_color_iter += COMPONENTS_PER_COLOR;
    }

    *unique_color_count = color_set.color_count;
    return unique_colors;
}

u8 *colors_unique_inplace(u8 *colors, isize color_count, isize *unique_color_count, void *arena) {
    if (color_count == 0) {
        *unique_color_count = 0;
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    ColorSet color_set;

    // Must be a positive integer in the for of (2^N - 1).
    isize const INITIAL_INDEX_COUNT = 63;
    i32 *initial_indices = arena_alloc(arena, INITIAL_INDEX_COUNT * sizeof(i32));
    if (initial_indices == NULL) {
        *unique_color_count = 0;
        return NULL;
    }
    memset(initial_indices, -1, INITIAL_INDEX_COUNT * sizeof(i32));

    color_set.index_count = INITIAL_INDEX_COUNT;
    for (int i = 0; i < USIZE_BITS - usize_leading_zeroes(INITIAL_INDEX_COUNT); i += 1) {
        color_set.indices[i] = &initial_indices[(1 << i) - 1];
    }

    color_set.color_count = 0;
    color_set.color_capacity = color_count;
    isize color_chunk, color_chunk_pos;
    chunk_pos(color_count - 1, &color_chunk, &color_chunk_pos);
    for (int i = 0; i <= color_chunk; i += 1) {
        color_set.colors[i] = &colors[((1 << i) - 1) * COMPONENTS_PER_COLOR];
    }

    if (!color_set_append(&color_set, colors, color_count, arena)) {
        *unique_color_count = 0;
        return NULL;
    }

    arena_rewind(arena, snapshot);
    *unique_color_count = color_set.color_count;
    return color_set.colors[0];
}

f32 *palette_by_median_cut(
    f32 const *colors_raw, isize color_count,
    isize target_color_count,
    isize *colors_generated,
    void *arena
) {
    ArenaSnapshot snapshot = arena_snapshot(arena);

    f32 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (palette == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    f32x3 *colors = arena_alloc(arena, color_count * sizeof(f32x3));
    if (colors == NULL) {
        *colors_generated = 0;
        return NULL;
    }
    memcpy(colors, colors_raw, color_count * sizeof(f32x3));

    typedef struct {
        f32x3 *begin;
        f32x3 *end;
    } Segment;

    // Ring buffer is used for the queue.
    typedef struct {
        Segment *segments;
        isize count;
        isize capacity;

        isize first_index;
        isize last_index;
    } SegmentsQueue;

    Segment *queue_segments = arena_alloc(arena, target_color_count * sizeof(Segment));
    if (queue_segments == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    SegmentsQueue queue = {
        .segments = queue_segments,
        .count = 1,
        .capacity = target_color_count,
        .first_index = 0,
        .last_index = 0,
    };
    queue.segments[queue.first_index] = (Segment){colors, colors + color_count};

    while (queue.count < target_color_count) {
        assert(queue.count > 0);
        Segment current_segment = queue.segments[queue.first_index];

        // The image probably has less colors then we are targeting. Break out of the loop because
        // segments are naturally sorted by their size in decreasing order, so the rest of the
        // segments are going to be of size 1 as well.
        if (current_segment.end - current_segment.begin == 1) {
            break;
        }

        queue.first_index = (queue.first_index + 1) % queue.capacity;
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
            comparator = f32x3_compare_by_x;
        } else if (max_y - min_y >= max_x - min_x && max_y - min_y >= max_z - min_z) {
            comparator = f32x3_compare_by_y;
        } else {
            comparator = f32x3_compare_by_z;
        }
        qsort(
            current_segment.begin, current_segment.end - current_segment.begin, sizeof(f32x3),
            comparator
        );

        queue.segments[(queue.last_index + 1) % queue.capacity] = (Segment){
            .begin = current_segment.begin,
            .end = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
        };
        queue.segments[(queue.last_index + 2) % queue.capacity] = (Segment){
            .begin = current_segment.begin + (current_segment.end - current_segment.begin) / 2,
            .end = current_segment.end,
        };
        queue.last_index = (queue.last_index + 2) % queue.capacity;
        queue.count += 2;
    }

    f32 *palette_iter = palette;
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
        palette_iter[0] = (f32)(x_sum / (f64)color_count);
        palette_iter[1] = (f32)(y_sum / (f64)color_count);
        palette_iter[2] = (f32)(z_sum / (f64)color_count);

        if (segment_index == queue.last_index) {
            break;
        }
        palette_iter += COMPONENTS_PER_COLOR;
        segment_index = (segment_index + 1) % queue.capacity;
    }

    // Tighten the initial palette allocation.
    arena_rewind(arena, snapshot);
    palette = arena_alloc(arena, queue.count * COMPONENTS_PER_COLOR * sizeof(f32));
    *colors_generated = queue.count;

    return palette;
}

f32 *palette_by_k_means(
    f32 const *colors_raw, isize color_count,
    isize target_color_count,
    isize *colors_generated,
    bool k_means_plus_plus,
    void *arena
) {
    isize const MAX_ITERATIONS = 300;

    ArenaSnapshot snapshot = arena_snapshot(arena);

    f32 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(f32));
    if (palette == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    f32x3 *colors = arena_alloc(arena, color_count * sizeof(f32x3));
    f32x3 *colors_end = colors + color_count;
    if (colors == NULL) {
        *colors_generated = 0;
        return NULL;
    }
    {
        f32 const *raw_color_iter = colors_raw;
        for (isize i = 0; i < color_count; i += 1) {
            colors[i] = (f32x3){raw_color_iter[0], raw_color_iter[1], raw_color_iter[2]};
            raw_color_iter += COMPONENTS_PER_COLOR;
        }
    }

    // In total the whole image contained less colors than we asked for.
    if (color_count <= target_color_count) {
        memcpy(palette, colors, color_count * sizeof(f32x3));
        *colors_generated = color_count;
        arena_rewind(arena, snapshot);
        return palette;
    }

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (u64)memcpy, (u64)memmove);

    // Pick the initial set of centroids randomly from the image colors.
    f32x3 *centroids = arena_alloc(arena, target_color_count * sizeof(f32x3));
    f32x3 *new_centroids = arena_alloc(arena, target_color_count * sizeof(f32x3));
    if (centroids == NULL || new_centroids == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    isize centroid_count = 0;
    if (k_means_plus_plus) {
        ArenaSnapshot snapshot = arena_snapshot(arena);
        f32 *centroid_distances = arena_alloc(arena, color_count * sizeof(f32));
        if (centroid_distances == NULL) {
            *colors_generated = 0;
            return NULL;
        }

        // Choose the first centroid purely randomly.
        {
            isize centroid_index = f64_random(&rng) * color_count;
            centroids[0] = colors[centroid_index];

            f32x3 swap = colors[centroid_index];
            colors[centroid_index] = colors[0];
            colors[0] = swap;

            centroid_count = 1;
        }

        // Sum of the distances to the closest centroids for the remaining colors.
        f64 centroid_distance_sum = 0.0;
        for (isize i = centroid_count; i < color_count; i += 1) {
             f32 distance = f32x3_dot(
                f32x3_sub(centroids[0], colors[i]),
                f32x3_sub(centroids[0], colors[i])
            );
            centroid_distances[i] = distance;
            centroid_distance_sum += distance;
        }

        while (centroid_count < target_color_count) {
            isize centroid_index;
            f32x3 centroid;
            {
                f64 random = f64_random(&rng) * centroid_distance_sum;
                f64 running_sum = 0.0;

                // Iterate from "centroid_count" index to skip colors already chosen as centroids.
                for (isize i = centroid_count; i < color_count; i += 1) {
                    running_sum += centroid_distances[i];
                    if (running_sum >= random) {
                        centroid_index = i;
                        centroid = colors[i];
                        break;
                    }
                }
            }

            centroids[centroid_count] = centroid;
            {
                f32 swap = centroid_distances[centroid_index];
                centroid_distances[centroid_index] = centroid_distances[centroid_count];
                centroid_distances[centroid_count] = swap;
            }
            {
                f32x3 swap = colors[centroid_index];
                colors[centroid_index] = colors[centroid_count];
                colors[centroid_count] = swap;
            }

            for (isize i = centroid_count; i < color_count; i += 1) {
                f32 distance = f32x3_dot(
                    f32x3_sub(colors[i], centroid),
                    f32x3_sub(colors[i], centroid)
                );
                if (distance < centroid_distances[i]) {
                    centroid_distance_sum -= centroid_distances[i] - distance;
                    centroid_distances[i] = distance;
                }
            }

            centroid_count += 1;
        }

        arena_rewind(arena, snapshot);
    } else {
        f32x3 *centroid_iter = centroids;
        f32x3 *color_iter = colors;
        for (isize i = 0; i < target_color_count; i += 1) {
            isize next_color_index = f64_random(&rng) * (colors_end - color_iter);
            assert(next_color_index < colors_end - color_iter);

            *centroid_iter = color_iter[next_color_index];
            centroid_iter += 1;

            f32x3 swap = color_iter[next_color_index];
            color_iter[next_color_index] = color_iter[0];
            color_iter[0] = swap;
            color_iter += 1;

            centroid_count += 1;
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

    // Tighten the initial palette allocation.
    arena_rewind(arena, snapshot);
    palette = arena_alloc(arena, centroid_count * COMPONENTS_PER_COLOR * sizeof(f32));
    memcpy(palette, centroids, centroid_count * sizeof(f32x3));
    *colors_generated = centroid_count;

    return palette;
}

u8 *palette_by_modified_median_cut(
    u8 const *colors, isize color_count,
    isize target_color_count,
    isize *colors_generated,
    void *arena
) {
    int const QUANT = 6;

    ArenaSnapshot snapshot = arena_snapshot(arena);

    u8 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (palette == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    // Quantize input colors to QUANT bits and count colors within each "quantum volume".
    // LSB of an index are for the red component, MSB are for the blue component.

    isize quantized_color_count = (isize)(1 << QUANT) * (1 << QUANT) * (1 << QUANT);
    isize *color_frequencies = arena_alloc(arena, quantized_color_count * sizeof(isize));
    if (color_frequencies == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    memset(color_frequencies, 0, quantized_color_count * sizeof(isize));
    {
        u8 const *color_iter = colors;
        while (color_iter < colors + color_count * COMPONENTS_PER_COLOR) {
            isize index =
                (u32)(color_iter[0] >> (8 - QUANT)) << (0 * QUANT) |
                (u32)(color_iter[1] >> (8 - QUANT)) << (1 * QUANT) |
                (u32)(color_iter[2] >> (8 - QUANT)) << (2 * QUANT);
            color_frequencies[index] += 1;

            color_iter += COMPONENTS_PER_COLOR;
        }
    }

    // The "min" end of the box is inclusive, the "max" end is exclusive.
    typedef struct {
        isize min[3];
        isize max[3];
        isize color_count;
        isize volume;
    } Box;

    typedef struct {
        Box *boxes;
        isize active_count;
        isize total_count;
        isize capacity;
    } BoxQueue;

    // We start with a single box encompassing the whole color space.
    // Presumably this way the generated palette works better with error diffusion dithering.
    BoxQueue queue = {
        .boxes = arena_alloc(arena, target_color_count * sizeof(Box)),
        .active_count = 1,
        .total_count = 1,
        .capacity = target_color_count,
    };
    if (queue.boxes == NULL) {
        *colors_generated = 0;
        return NULL;
    }
    queue.boxes[0] = (Box){
        .min = {0, 0, 0},
        .max = {1 << QUANT, 1 << QUANT, 1 << QUANT},
        .color_count = color_count,
        .volume = (isize)(1 << QUANT) * (1 << QUANT) * (1 << QUANT),
    };

    isize *box_plane_frequencies = arena_alloc(arena, (1 << QUANT) * sizeof(isize));
    if (box_plane_frequencies == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    // First choose the next box to subdivide based on color count, do the rest of the subdivisions
    // based on combination of both volume and color count.
    bool volume_affects_priority = false;

    while (queue.active_count > 0 && queue.total_count < target_color_count) {
        // Once we've generated 80% of the colors, start counting box volume into account as well.
        if (
            !volume_affects_priority &&
            (target_color_count - queue.total_count) * 5 < target_color_count * 4
        ) {
            volume_affects_priority = true;

            // Re-heapify the queue based on the changed priority metric.
            for (isize i = 1; i < queue.active_count; i += 1) {
                isize current = i;
                Box current_box = queue.boxes[current];
                isize current_priority = current_box.color_count * current_box.volume;

                while (current > 0) {
                    isize parent = (current - 1) / 2;
                    isize parent_priority =
                        queue.boxes[parent].color_count * queue.boxes[parent].volume;

                    if (parent_priority > current_priority) {
                        break;
                    }
                    queue.boxes[current] = queue.boxes[parent];
                    current = parent;
                }

                queue.boxes[current] = current_box;
            }
        }

        Box next_box = queue.boxes[0];

        // Pick the longest dimension of the box.
        int this_dimension = 0;
        if (
            next_box.max[1] - next_box.min[1] >= next_box.max[0] - next_box.min[0] &&
            next_box.max[1] - next_box.min[1] >= next_box.max[2] - next_box.min[2]
        ) {
            this_dimension = 1;
        } else if (
            next_box.max[2] - next_box.min[2] >= next_box.max[0] - next_box.min[0] &&
            next_box.max[2] - next_box.min[2] >= next_box.max[1] - next_box.min[1]
        ) {
            this_dimension = 2;
        }
        int other_dimension = (this_dimension + 1) % 3;
        int another_dimension = (this_dimension + 2) % 3;

        // Find the median point.
        for (
            isize this = next_box.min[this_dimension];
            this < next_box.max[this_dimension];
            this += 1
        ) {
            box_plane_frequencies[this] = 0;

            for (
                isize other = next_box.min[other_dimension];
                other < next_box.max[other_dimension];
                other += 1
            ) {
                for (
                    isize another = next_box.min[another_dimension];
                    another < next_box.max[another_dimension];
                    another += 1
                ) {
                    isize index =
                        this    << (this_dimension    * QUANT) |
                        other   << (other_dimension   * QUANT) |
                        another << (another_dimension * QUANT);

                    box_plane_frequencies[this] += color_frequencies[index];
                }
            }
        }
        isize median_index = next_box.min[this_dimension];
        isize left_box_color_count = 0;
        while (
            median_index < next_box.max[this_dimension] &&
            left_box_color_count <= next_box.color_count / 2
        ) {
            left_box_color_count += box_plane_frequencies[median_index];
            median_index += 1;
        }
        {
            // Put the box at the median point into the smaller half.
            isize left_half = median_index - next_box.min[this_dimension];
            isize right_half = next_box.max[this_dimension] - median_index;
            if (left_half < right_half) {
                left_box_color_count += box_plane_frequencies[median_index];
                median_index += 1;
            }
        }

        // Move the median partition point further into the less dense partition.
        isize partition_index;
        {
            isize left_half = median_index - next_box.min[this_dimension];
            isize right_half = next_box.max[this_dimension] - median_index;
            if (left_half > right_half) {
                partition_index = isize_clamp(
                    (next_box.min[this_dimension] + median_index) / 2,
                    next_box.min[this_dimension] + 1,
                    next_box.max[this_dimension] - 1
                );
                for (isize i = partition_index; i < median_index; i += 1) {
                    left_box_color_count -= box_plane_frequencies[i];
                }
            } else if (left_half < right_half) {
                partition_index = isize_clamp(
                    (median_index + next_box.max[this_dimension]) / 2,
                    next_box.min[this_dimension] + 1,
                    next_box.max[this_dimension] - 1
                );
                for (isize i = median_index; i < partition_index; i += 1) {
                    left_box_color_count += box_plane_frequencies[i];
                }
            } else {
                partition_index = isize_clamp(
                    median_index,
                    next_box.min[this_dimension] + 1,
                    next_box.max[this_dimension] - 1
                );
            }
        }

        // Split the "next_box" into left and right boxes.

        Box left_box = next_box;
        left_box.max[this_dimension] = partition_index;
        left_box.color_count = left_box_color_count;
        left_box.volume =
            (left_box.max[0] - left_box.min[0]) *
            (left_box.max[1] - left_box.min[1]) *
            (left_box.max[2] - left_box.min[2]);

        Box right_box = next_box;
        right_box.min[this_dimension] = partition_index;
        right_box.color_count = next_box.color_count - left_box.color_count;
        right_box.volume = next_box.volume - left_box.volume;

        // We've reached the quantum volume box.
        // Can't subdivide it any further, put it away at the very and of the queue buffer.
        if (next_box.max[this_dimension] - next_box.min[this_dimension] == 1) {
            queue.boxes[queue.capacity - (queue.total_count - queue.active_count) - 1] = next_box;
            queue.active_count -= 1;
            left_box = queue.boxes[queue.active_count];
            right_box = (Box){0};
        }

        // Replace the "next_box" with the left box in the queue.
        if (left_box.max[this_dimension] - left_box.min[this_dimension] > 0) {
            isize current = 0;
            isize current_priority = volume_affects_priority
                ? left_box.color_count * left_box.volume
                : left_box.color_count;

            while (2 * current + 1 < queue.active_count) {
                isize left = 2 * current + 1;
                isize left_priority = volume_affects_priority
                    ? queue.boxes[left].color_count * queue.boxes[left].volume
                    : queue.boxes[left].color_count;

                isize right = 2 * current + 2;
                isize right_priority = -1;
                if (right < queue.active_count) {
                    right_priority = volume_affects_priority
                        ? queue.boxes[right].color_count * queue.boxes[right].volume
                        : queue.boxes[right].color_count;
                }

                isize max = current;
                isize max_priority = current_priority;

                if (left_priority > max_priority) {
                    max = left;
                    max_priority = left_priority;
                }
                if (right_priority > max_priority) {
                    max = right;
                    max_priority = right_priority;
                }

                if (current == max) {
                    break;
                }
                queue.boxes[current] = queue.boxes[max];
                current = max;
            }

            queue.boxes[current] = left_box;
        }

        // Add the right box to the queue.
        if (right_box.max[this_dimension] - right_box.min[this_dimension] > 0) {
            queue.active_count += 1;
            queue.total_count += 1;

            isize current = queue.active_count - 1;
            isize current_priority = volume_affects_priority
                ? right_box.color_count * right_box.volume
                : right_box.color_count;

            while (current > 0) {
                isize parent = (current - 1) / 2;
                isize parent_priority = volume_affects_priority
                    ? queue.boxes[parent].color_count * queue.boxes[parent].volume
                    : queue.boxes[parent].color_count;

                if (parent_priority > current_priority) {
                    break;
                }
                queue.boxes[current] = queue.boxes[parent];
                current = parent;
            }

            queue.boxes[current] = right_box;
        }
    }

    u8 *palette_iter = palette;
    for (isize i = 0; i < queue.total_count; i += 1) {
        Box box;
        if (i < queue.active_count) {
            box = queue.boxes[i];
        } else {
            // Coming back for those quantum volume boxes which we couldn't subdivide any further.
            box = queue.boxes[queue.capacity - i - 1];
        }

        f64 sum[3] = {0};
        if (box.color_count > 0) {
            for (isize i = box.min[0]; i < box.max[0]; i += 1) {
                for (isize j = box.min[1]; j < box.max[1]; j += 1) {
                    for (isize k = box.min[2]; k < box.max[2]; k += 1) {
                        isize index = i | j << QUANT | k << QUANT << QUANT;
                        sum[0] += (i << (8 - QUANT)) / 255.0F * color_frequencies[index];
                        sum[1] += (j << (8 - QUANT)) / 255.0F * color_frequencies[index];
                        sum[2] += (k << (8 - QUANT)) / 255.0F * color_frequencies[index];
                    }
                }
            }
        }

        if (box.color_count == 0) {
            palette_iter[0] = ((box.min[0] << (8 - QUANT)) + ((box.max[0] - 1) << (8 - QUANT))) / 2;
            palette_iter[1] = ((box.min[1] << (8 - QUANT)) + ((box.max[1] - 1) << (8 - QUANT))) / 2;
            palette_iter[2] = ((box.min[2] << (8 - QUANT)) + ((box.max[2] - 1) << (8 - QUANT))) / 2;
        } else {
            palette_iter[0] = sum[0] / box.color_count * 255.0F;
            palette_iter[1] = sum[1] / box.color_count * 255.0F;
            palette_iter[2] = sum[2] / box.color_count * 255.0F;
        }

        palette_iter += COMPONENTS_PER_COLOR;
    }

    // Tighten the initial palette allocation.
    arena_rewind(arena, snapshot);
    palette = arena_alloc(arena, queue.total_count * COMPONENTS_PER_COLOR * sizeof(u8));
    *colors_generated = queue.total_count;

    return palette;
}

u8 *palette_by_octree(
    u8 const *colors, isize color_count,
    isize target_color_count,
    isize *colors_generated,
    void *arena
) {
    ArenaSnapshot snapshot = arena_snapshot(arena);

    u8 *palette = arena_alloc(arena, target_color_count * COMPONENTS_PER_COLOR * sizeof(u8));
    if (palette == NULL) {
        *colors_generated = 0;
        return NULL;
    }

    enum { MAX_LEVEL = 8 };

    typedef struct Node Node;

    struct Node {
        u64 color[3];
        isize reference_count;
        Node *children[8];

        Node *next_sibling;
        Node *previous_sibling;
    };

    Node root = {
        .color = {0},
        .reference_count = 0,
    };

    struct {
        Node *first;
        Node *last;
    } nodes_by_level[MAX_LEVEL + 1] = {[0] = {&root, &root}};

    isize leaf_node_count = 0;

    // Push the colors to the very bottom of the octree.
    u8 const *color_iter = colors;
    while (color_iter < colors + color_count * COMPONENTS_PER_COLOR) {
        Node *tree_iter = &root;
        isize level = 0;
        while (level < MAX_LEVEL) {
            isize child_index =
                ((color_iter[0] >> (MAX_LEVEL - level - 1)) & 0x1) << 0 |
                ((color_iter[1] >> (MAX_LEVEL - level - 1)) & 0x1) << 1 |
                ((color_iter[2] >> (MAX_LEVEL - level - 1)) & 0x1) << 2;

            if (tree_iter->children[child_index] == NULL) {
                Node *new_child = arena_alloc(arena, sizeof(Node));
                if (new_child == NULL) {
                    *colors_generated = 0;
                    return NULL;
                }
                memset(new_child, 0, sizeof(Node));

                // "level + 1" is because the current node children are on the level below.

                if (nodes_by_level[level + 1].first == NULL) {
                    nodes_by_level[level + 1].first = new_child;
                    nodes_by_level[level + 1].last = new_child;
                } else {
                    nodes_by_level[level + 1].last->next_sibling = new_child;
                    new_child->previous_sibling = nodes_by_level[level + 1].last;

                    nodes_by_level[level + 1].last = new_child;
                }

                if (level + 1 == MAX_LEVEL) {
                    leaf_node_count += 1;
                }

                tree_iter->children[child_index] = new_child;
            }

            tree_iter = tree_iter->children[child_index];
            level += 1;
        }

        tree_iter->reference_count += 1;
        tree_iter->color[0] += color_iter[0];
        tree_iter->color[1] += color_iter[1];
        tree_iter->color[2] += color_iter[2];

        color_iter += COMPONENTS_PER_COLOR;
    }

    // Reduce the tree nodes until we reach the required number of leaves in the tree.
    // Start merging from one level above the very bottom one.
    isize level = MAX_LEVEL;
    Node *sibling_iter = NULL;
    while (leaf_node_count > target_color_count) {
        if (sibling_iter == NULL) {
            level -= 1;
            sibling_iter = nodes_by_level[level].first;
        }
        if (level < 0) {
            break;
        }

        for (isize i = 0; i < countof(sibling_iter->children); i += 1) {
            Node *child = sibling_iter->children[i];

            if (child != NULL) {
                sibling_iter->reference_count += child->reference_count;
                sibling_iter->color[0] += child->color[0];
                sibling_iter->color[1] += child->color[1];
                sibling_iter->color[2] += child->color[2];

                if (child->previous_sibling != NULL) {
                    child->previous_sibling->next_sibling = child->next_sibling;
                }
                if (child->next_sibling != NULL) {
                    child->next_sibling->previous_sibling = child->previous_sibling;
                }
                if (child == nodes_by_level[level + 1].first) {
                    nodes_by_level[level + 1].first = child->next_sibling;
                }
                if (child == nodes_by_level[level + 1].last) {
                    nodes_by_level[level + 1].last = child->previous_sibling;
                }

                leaf_node_count -= 1;
            }

            sibling_iter->children[i] = NULL;
        }

        // The node that we merged the children into became the leaf node itself.
        leaf_node_count += 1;

        sibling_iter = sibling_iter->next_sibling;
    }

    u8 *palette_iter = palette;
    isize palette_size = 0;

    level = MAX_LEVEL + 1;
    sibling_iter = NULL;
    for (isize i = 0; i < target_color_count; i += 1) {
        while (sibling_iter == NULL && level >= 0) {
            level -= 1;
            sibling_iter = nodes_by_level[level].first;
        }
        if (level < 0) {
            break;
        }

        if (sibling_iter->reference_count == 0) {
            break;
        }

        palette_iter[0] = sibling_iter->color[0] / sibling_iter->reference_count;
        palette_iter[1] = sibling_iter->color[1] / sibling_iter->reference_count;
        palette_iter[2] = sibling_iter->color[2] / sibling_iter->reference_count;

        sibling_iter = sibling_iter->next_sibling;

        palette_iter += COMPONENTS_PER_COLOR;
        palette_size += 1;
    }

    // Tighten the initial palette allocation.
    arena_rewind(arena, snapshot);
    *colors_generated = palette_size;
    palette = arena_alloc(arena, palette_size * COMPONENTS_PER_COLOR * sizeof(f32));

    return palette;
}

#define COLOR_TABLE_GRID_SIZE 16

#define COLOR_TABLE_CELL_COUNT  \
    ((isize)                    \
        COLOR_TABLE_GRID_SIZE * \
        COLOR_TABLE_GRID_SIZE * \
        COLOR_TABLE_GRID_SIZE   \
    )

typedef struct {
    f32x3 color;
    GifColorIndex index;
} IndexedColor;

// A list of closest-color candidates for a grid cell to brute-force on.
typedef struct {
    isize count;
    IndexedColor colors[];
} ColorTableCandidates;

typedef struct {
    f32xbox bounding_box;
    f32x3 cell_size;
    ColorTableCandidates **candidates;
} ColorTable;

static ColorTable *color_table_create(f32 const *colors, isize color_count, Arena *arena) {
    assert(0 < color_count && color_count <= (isize)GIF_COLOR_INDEX_MAX + 1);

    ColorTable *table = arena_alloc(arena, sizeof(ColorTable));
    if (table == NULL) {
        return NULL;
    }

    // A bounding box for the "colors" array.
    {
        table->bounding_box.min = (f32x3){ INFINITY,  INFINITY,  INFINITY};
        table->bounding_box.max = (f32x3){-INFINITY, -INFINITY, -INFINITY};

        f32 const *color_iter = colors;
        while (color_iter < colors + COMPONENTS_PER_COLOR * color_count) {
            f32x3 color = {color_iter[0], color_iter[1], color_iter[2]};
            table->bounding_box.min = f32x3_min(table->bounding_box.min, color);
            table->bounding_box.max = f32x3_max(table->bounding_box.max, color);

            color_iter += COMPONENTS_PER_COLOR;
        }
    }
    table->cell_size = f32x3_scale(
        f32x3_sub(table->bounding_box.max, table->bounding_box.min),
        1.0F / COLOR_TABLE_GRID_SIZE
    );

    table->candidates = arena_alloc(arena, COLOR_TABLE_CELL_COUNT * sizeof(ColorTableCandidates *));
    if (table->candidates == NULL) {
        return NULL;
    }

    // Create a list of candidate palette colors for each cell. When an image pixel falls into some
    // cell, the closest palette color can be found among the candidates.
    for (isize cell_index = 0; cell_index < COLOR_TABLE_CELL_COUNT; cell_index += 1) {
        f32x3 cell_coords = {
            cell_index
                % COLOR_TABLE_GRID_SIZE,
            cell_index / COLOR_TABLE_GRID_SIZE
                % COLOR_TABLE_GRID_SIZE,
            cell_index / (COLOR_TABLE_GRID_SIZE * COLOR_TABLE_GRID_SIZE)
                % COLOR_TABLE_GRID_SIZE,
        };

        f32xbox cell;
        cell.min = f32x3_add(table->bounding_box.min, f32x3_mul(table->cell_size, cell_coords));
        cell.max = f32x3_add(cell.min, table->cell_size);

        // Find the "best" (i.e. min) distance out of all possible max distances from the palette
        // colors to the points inside the cell. If a min distance from some color to the cell
        // points is larger than that, then this color is not worth checking.
        f32 best_max_distance = INFINITY;
        {
            f32 const *color_iter = colors;
            while (color_iter < colors + color_count * COMPONENTS_PER_COLOR) {
                f32x3 color = {color_iter[0], color_iter[1], color_iter[2]};

                // Obviously works when the color component is between min and max.
                // Otherwise works because one of the differences gets negative.
                f32x3 farthest_point = {
                    color.x - cell.min.x >= cell.max.x - color.x ? cell.min.x : cell.max.x,
                    color.y - cell.min.y >= cell.max.y - color.y ? cell.min.y : cell.max.y,
                    color.z - cell.min.z >= cell.max.z - color.z ? cell.min.z : cell.max.z,
                };
                f32 max_distance = f32x3_dot(
                    f32x3_sub(color, farthest_point),
                    f32x3_sub(color, farthest_point)
                );

                if (max_distance < best_max_distance) {
                    best_max_distance = max_distance;
                }

                color_iter += COMPONENTS_PER_COLOR;
            }
        }

        // Conservative allocation in front, tighten it later. Allocating memory this way also plays
        // nicely with an idea of extending the arena memory using the GIF_LIB_REPORT_ALLOC hook.
        ArenaSnapshot snapshot = arena_snapshot(arena);

        ColorTableCandidates *candidates = arena_alloc(
            arena,
            sizeof(ColorTableCandidates) + color_count * sizeof(IndexedColor)
        );
        if (candidates == NULL) {
            return NULL;
        }
        candidates->count = 0;

        // Populate the list of candidate colors with the ones worth checking when a query color
        // falls into the current cell.
        {
            GifColorIndex color_index = 0;
            f32 const *color_iter = colors;

            while (color_iter < colors + color_count * COMPONENTS_PER_COLOR) {
                f32x3 color = {color_iter[0], color_iter[1], color_iter[2]};

                f32x3 closest_point = f32x3_clamp(color, cell.min, cell.max);
                f32 min_distance = f32x3_dot(
                    f32x3_sub(color, closest_point),
                    f32x3_sub(color, closest_point)
                );

                // Otherwise this color is not worth looking into: a better candidate for any point
                // inside the current cell is the one which had "best_max_distance" distance in the
                // worst case.
                if (min_distance <= best_max_distance) {
                    candidates->colors[candidates->count++] = (IndexedColor){
                        .color = color,
                        .index = color_index,
                    };
                }

                color_iter += COMPONENTS_PER_COLOR;
                color_index += 1;
            }
        }

        arena_rewind(arena, snapshot);
        table->candidates[cell_index] = arena_alloc(
            arena,
            sizeof(ColorTableCandidates) + candidates->count * sizeof(IndexedColor)
        );
    }

    return table;
}

static IndexedColor color_table_get_closest(ColorTable const *table, f32x3 search_color) {
    // If a color falls outside of the palette bounding box, look through the candidates of a
    // clamped point. I'm not exactly sure, but I think there is no way a closer color could be
    // found among the candidates from the neighboring cells (or any other cells).
    isize cell_x = isize_clamp(
        (search_color.x - table->bounding_box.min.x) / table->cell_size.x,
        0, COLOR_TABLE_GRID_SIZE - 1
    );
    isize cell_y = isize_clamp(
        (search_color.y - table->bounding_box.min.y) / table->cell_size.y,
        0, COLOR_TABLE_GRID_SIZE - 1
    );
    isize cell_z = isize_clamp(
        (search_color.z - table->bounding_box.min.z) / table->cell_size.z,
        0, COLOR_TABLE_GRID_SIZE - 1
    );

    isize cell_index =
        cell_x +
        cell_y * COLOR_TABLE_GRID_SIZE +
        cell_z * COLOR_TABLE_GRID_SIZE * COLOR_TABLE_GRID_SIZE;

    IndexedColor nearest;
    f32 distance_to_nearest = INFINITY;
    ColorTableCandidates *candidates = table->candidates[cell_index];

    for (isize i = 0; i < candidates->count; i += 1) {
        IndexedColor candidate = candidates->colors[i];
        f32 distance = f32x3_dot(
            f32x3_sub(candidate.color, search_color),
            f32x3_sub(candidate.color, search_color)
        );

        if (distance < distance_to_nearest) {
            distance_to_nearest = distance;
            nearest = candidate;
        }
    }

    return nearest;
}

GifColorIndex *image_quantize(
    f32 const *pixels, isize pixel_count,
    f32 const *colors, isize color_count,
    void *arena
) {
    assert(0 < color_count && color_count <= (isize)GIF_COLOR_INDEX_MAX + 1);

    GifColorIndex *indexed_pixels = arena_alloc(arena, pixel_count * sizeof(GifColorIndex));
    if (indexed_pixels == NULL) {
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    ColorTable *table = color_table_create(colors, color_count, arena);
    if (table == NULL) {
        return NULL;
    }

    f32 const *pixel_iter = pixels;
    GifColorIndex *indexed_pixel_iter = indexed_pixels;
    while (pixel_iter != pixels + pixel_count * COMPONENTS_PER_COLOR) {
        f32x3 pixel = {pixel_iter[0], pixel_iter[1], pixel_iter[2]};
        *indexed_pixel_iter = color_table_get_closest(table, pixel).index;

        pixel_iter += COMPONENTS_PER_COLOR;
        indexed_pixel_iter += 1;
    }

    arena_rewind(arena, snapshot);
    return indexed_pixels;
}

GifColorIndex *image_floyd_steinberg_dither(
    f32 const *image, isize width, isize height,
    f32 const *colors, isize color_count,
    void *arena
) {
    assert(0 < color_count && color_count <= (isize)GIF_COLOR_INDEX_MAX + 1);

    GifColorIndex *indexed_pixels = arena_alloc(arena, width * height * sizeof(GifColorIndex));
    if (indexed_pixels == NULL) {
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    // Copy the input image into a separate array, because we will be modifying the input pixels.
    f32x3 *pixels = arena_alloc(arena, width * height * sizeof(f32x3));
    if (pixels == NULL) {
        return NULL;
    }
    {
        f32 const *image_iter = image;
        for (isize i = 0; i < width * height; i += 1) {
            pixels[i] = (f32x3){image_iter[0], image_iter[1], image_iter[2]};
            image_iter += COMPONENTS_PER_COLOR;
        }
    }

    ColorTable *table = color_table_create(colors, color_count, arena);
    if (table == NULL) {
        return NULL;
    }

    // Limit error diffusion with an arbitrary value of 5% of the color palette bounding box span.
    f32xbox const error_clamp = {
        .min = f32x3_scale(f32x3_sub(table->bounding_box.min, table->bounding_box.max), 0.05F),
        .max = f32x3_scale(f32x3_sub(table->bounding_box.max, table->bounding_box.min), 0.05F),
    };

    // Limit the pixel values. I'm not sure if it's even necessary to clamp the pixel values.
    f32xbox const pixel_clamp = {
        .min = f32x3_add(table->bounding_box.min, error_clamp.min),
        .max = f32x3_add(table->bounding_box.max, error_clamp.max),
    };

    for (isize y = 0; y < height; y += 1) {
        for (isize x = 0; x < width; x += 1) {
            f32x3 current_color = pixels[y * width + x];
            IndexedColor closest = color_table_get_closest(table, current_color);

            indexed_pixels[y * width + x] = closest.index;

            f32x3 error = f32x3_clamp(
                f32x3_sub(current_color, closest.color),
                error_clamp.min, error_clamp.max
            );

            if (x + 1 < width) {
                pixels[y * width + (x + 1)] = f32x3_clamp(
                    f32x3_add(pixels[y * width + (x + 1)], f32x3_scale(error, 7.0F / 16.0F)),
                    pixel_clamp.min, pixel_clamp.max
                );
            }

            if (x - 1 >= 0 && y + 1 < height) {
                pixels[(y + 1) * width + (x - 1)] = f32x3_clamp(
                    f32x3_add(pixels[(y + 1) * width + (x - 1)], f32x3_scale(error, 3.0F / 16.0F)),
                    pixel_clamp.min, pixel_clamp.max
                );
            }

            if (y + 1 < height) {
                pixels[(y + 1) * width + x] = f32x3_clamp(
                    f32x3_add(pixels[(y + 1) * width + x], f32x3_scale(error, 5.0F / 16.0F)),
                    pixel_clamp.min, pixel_clamp.max
                );
            }

            if (x + 1 < width && y + 1 < height) {
                pixels[(y + 1) * width + (x + 1)] = f32x3_clamp(
                    f32x3_add(pixels[(y + 1) * width + (x + 1)], f32x3_scale(error, 1.0F / 16.0F)),
                    pixel_clamp.min, pixel_clamp.max
                );
            }
        }
    }

    arena_rewind(arena, snapshot);
    return indexed_pixels;
}

GifColorIndex *image_ordered_dither(
    f32 const *image, isize width, isize height,
    f32 const *colors, isize color_count,
    void *arena
) {
    assert(0 < color_count && color_count <= (isize)GIF_COLOR_INDEX_MAX + 1);

    GifColorIndex *indexed_pixels = arena_alloc(arena, width * height * sizeof(GifColorIndex));
    if (indexed_pixels == NULL) {
        return NULL;
    }

    ArenaSnapshot snapshot = arena_snapshot(arena);

    ColorTable *table = color_table_create(colors, color_count, arena);
    if (table == NULL) {
        return NULL;
    }

    // https://en.wikipedia.org/wiki/Ordered_dithering
    // Generate the Bayer matrix.

    isize const BAYER_MATRIX_SIZE = 32;

    isize bayer_matrix_initial[4] = {
        0, 2,
        3, 1,
    };

    isize *bayer_matrix = bayer_matrix_initial;
    isize bayer_matrix_size = 2;

    while (bayer_matrix_size < BAYER_MATRIX_SIZE) {
        isize next_size = bayer_matrix_size * 2;
        isize *next_matrix = arena_alloc(arena, next_size * next_size * sizeof(isize));
        if (next_matrix == NULL) {
            return NULL;
        }

        for (isize y = 0; y < bayer_matrix_size; y += 1) {
            for (isize x = 0; x < bayer_matrix_size; x += 1) {
                isize previous = 4 * bayer_matrix[y * bayer_matrix_size + x];

                next_matrix[(y                ) * next_size + (x                )] = previous;
                next_matrix[(y                ) * next_size + (x + next_size / 2)] = previous + 2;
                next_matrix[(y + next_size / 2) * next_size + (x                )] = previous + 3;
                next_matrix[(y + next_size / 2) * next_size + (x + next_size / 2)] = previous + 1;
            }
        }

        bayer_matrix = next_matrix;
        bayer_matrix_size = next_size;
    }

    // Calculate the maximum distance between successive values on each channel in the palette.

    f32x3 *sorted_colors = arena_alloc(arena, color_count * sizeof(f32x3));
    if (sorted_colors == NULL) {
        return NULL;
    }
    {
        f32 const *color_iter = colors;
        for (isize i = 0; i < color_count; i += 1) {
            sorted_colors[i] = (f32x3){color_iter[0], color_iter[1], color_iter[2]};
            color_iter += COMPONENTS_PER_COLOR;
        }
    }

    f32x3 color_thresholds = {0.0F, 0.0F, 0.0F};

    qsort(sorted_colors, color_count, sizeof(f32x3), f32x3_compare_by_x);
    for (isize i = 1; i < color_count; i += 1) {
        color_thresholds.x = f32_max(
            color_thresholds.x,
            sorted_colors[i].x - sorted_colors[i - 1].x
        );
    }

    qsort(sorted_colors, color_count, sizeof(f32x3), f32x3_compare_by_y);
    for (isize i = 1; i < color_count; i += 1) {
        color_thresholds.y = f32_max(
            color_thresholds.y,
            sorted_colors[i].y - sorted_colors[i - 1].y
        );
    }

    qsort(sorted_colors, color_count, sizeof(f32x3), f32x3_compare_by_z);
    for (isize i = 1; i < color_count; i += 1) {
        color_thresholds.z = f32_max(
            color_thresholds.z,
            sorted_colors[i].z - sorted_colors[i - 1].z
        );
    }

    // Calculate a matrix with color offsets for each component.

    isize const OFFSET_MATRIX_SIZE = BAYER_MATRIX_SIZE;
    f32x3 *offset_matrix = arena_alloc(
        arena,
        OFFSET_MATRIX_SIZE * OFFSET_MATRIX_SIZE * sizeof(f32x3)
    );
    if (offset_matrix == NULL) {
        return NULL;
    }
    for (isize i = 0; i < OFFSET_MATRIX_SIZE * OFFSET_MATRIX_SIZE; i += 1) {
        f32 threshold = (f32)bayer_matrix[i] / (bayer_matrix_size * bayer_matrix_size);

        offset_matrix[i].x = color_thresholds.x * (threshold - 0.5F);
        offset_matrix[i].y = color_thresholds.y * (threshold - 0.5F);
        offset_matrix[i].z = color_thresholds.z * (threshold - 0.5F);
    }

    for (isize y = 0; y < height; y += 1) {
        for (isize x = 0; x < width; x += 1) {
            f32x3 current_color = {
                image[y * (width * COMPONENTS_PER_COLOR) + x * COMPONENTS_PER_COLOR + 0],
                image[y * (width * COMPONENTS_PER_COLOR) + x * COMPONENTS_PER_COLOR + 1],
                image[y * (width * COMPONENTS_PER_COLOR) + x * COMPONENTS_PER_COLOR + 2],
            };
            f32x3 shifted_color = f32x3_add(
                current_color,
                offset_matrix[y % OFFSET_MATRIX_SIZE * OFFSET_MATRIX_SIZE + x % OFFSET_MATRIX_SIZE]
            );

            indexed_pixels[y * width + x] = color_table_get_closest(table, shifted_color).index;
        }
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
    memset(out_buffer->data, 0, capacity);

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
    memset(new_buffer + out_buffer->capacity, 0, new_capacity - out_buffer->capacity);

    out_buffer->data = new_buffer;
    out_buffer->capacity = new_capacity;

    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    return true;
}

void gif_out_buffer_reset(GifOutputBuffer *out_buffer) {
    isize unencoded_size =
        (out_buffer->byte_pos - out_buffer->encoded_size) +
        (out_buffer->bit_pos == 0 ? 0 : 1);
    memmove(out_buffer->data, out_buffer->data + out_buffer->encoded_size, unencoded_size);

    out_buffer->byte_pos -= out_buffer->encoded_size;
    out_buffer->encoded_size = 0;

    u8 *free_memory_begin =
        out_buffer->data + out_buffer->byte_pos + (out_buffer->bit_pos == 0 ? 0 : 1);
    u8 *free_memory_end = out_buffer->data + out_buffer->capacity;
    memset(free_memory_begin, 0, free_memory_end - free_memory_begin);
}

typedef isize LzwCode;
#define LZW_CODE_INVALID ((LzwCode)(-1))
#define MAX_LZW_CODE ((LzwCode)4095)
#define LZW_CODE_MAX_BIT_LENGTH 12

typedef struct LzwTreeNode LzwTreeNode;

// A prefix tree for storing the "Color Sequence -> LZW Code" mappings.
// The first tree level (the root node itself) is for storing LZW codes for sequences of length one.
// LZW_CODE_INVALID is used to indicate an empty link.
// "(lzw_codes[i] == LZW_CODE_INVALID) => (children[i] == NULL)" is an invariant.
struct LzwTreeNode {
    LzwCode lzw_codes[GIF_MAX_COLORS];
    LzwTreeNode *children[GIF_MAX_COLORS];

    // Used for collecting tree nodes into "free" and "used" lists. (See "free_lzw_nodes" below.)
    LzwTreeNode *next;
};

typedef enum {
    GIF_ENCODER_IDLE,
    GIF_ENCODER_READY_FOR_NEXT_FRAME,
    GIF_ENCODER_FRAME_STARTED,
} GifEncoderState;

struct GifEncoder {
    GifEncoderState state;
    Arena *arena;

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

    struct {
        LzwTreeNode *root;
        LzwTreeNode *iter;

        // Last node added into the tree. The root is always the first one.
        LzwTreeNode *last;
    } lzw_tree;
    LzwCode next_lzw_code;

    // Nodes recycled after resetting the LZW tree.
    LzwTreeNode *free_lzw_nodes;
};

static bool gif_encoder_init(GifEncoder *encoder, Arena *arena) {
    encoder->state = GIF_ENCODER_IDLE;
    encoder->arena = arena;
    encoder->global_color_count = 0;
    encoder->local_color_count = 0;

    LzwTreeNode *root = arena_alloc(arena, sizeof(LzwTreeNode));
    if (root == NULL) {
        return false;
    }
    memset(root->children, 0, sizeof(root->children));

    // Initially fill the tree with every single color sequence of length one.
    for (isize index = 0; index < countof(root->lzw_codes); index += 1) {
        root->lzw_codes[index] = index;
    }

    encoder->free_lzw_nodes = NULL;
    encoder->lzw_tree.root = root;
    encoder->lzw_tree.last = root;

    // This means that we haven't started matching any color sequence from the LZW tree yet.
    encoder->lzw_tree.iter = root;
    encoder->next_lzw_code = LZW_CODE_INVALID;

    return true;
}

GifEncoder *gif_encoder_create(void *arena) {
    GifEncoder *encoder = arena_alloc(arena, sizeof(GifEncoder));
    if (encoder == NULL) {
        return NULL;
    }

    if (!gif_encoder_init(encoder, arena)) {
        return NULL;
    }

    // Amounts to a delay of 1 second.
    encoder->frame_delay = 100;

    return encoder;
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

        // Accumulate pixels from the image until we find a sequence which is not yet in the table.
        bool new_sequence_found;
        while (pixel_iter < pixels_end) {
            assert(*pixel_iter < color_count);

            new_sequence_found = encoder->lzw_tree.iter->lzw_codes[*pixel_iter] < 0;
            if (new_sequence_found) {
                // Don't consume the pixel which made our sequence "new".
                // On the next iteration the next accumulated sequence will start with this pixel.
                break;
            }

            if (encoder->lzw_tree.iter->children[*pixel_iter] == NULL) {
                LzwTreeNode *new_node;
                if (encoder->free_lzw_nodes != NULL) {
                    new_node = encoder->free_lzw_nodes;
                    encoder->free_lzw_nodes = encoder->free_lzw_nodes->next;
                } else {
                    new_node = arena_alloc(encoder->arena, sizeof(LzwTreeNode));
                    if (new_node == NULL) {
                        return -1;
                    }
                }
                encoder->lzw_tree.last->next = new_node;
                encoder->lzw_tree.last = new_node;
                memset(new_node->children, 0, sizeof(new_node->children));
                memset(new_node->lzw_codes, -1, sizeof(new_node->lzw_codes));

                encoder->lzw_tree.iter->children[*pixel_iter] = new_node;
            }

            encoder->next_lzw_code = encoder->lzw_tree.iter->lzw_codes[*pixel_iter];
            encoder->lzw_tree.iter = encoder->lzw_tree.iter->children[*pixel_iter];

            pixel_iter += 1;
        }

        if (new_sequence_found) {
            gif_encoder_write_lzw_code(
                encoder,
                encoder->next_lzw_code,
                encoder->lzw_code_bit_length,
                out_buffer
            );

            // Increase the LZW code length here because the spec says so:
            // > Whenever the LZW code value would exceed the current code length,
            // > the code length is increased by one. The packing/unpacking of these
            // > codes must then be altered to reflect the new code length.

            if (isize_power_of_two(encoder->next_available_lzw_code)) {
                encoder->lzw_code_bit_length += 1;
            }
            encoder->lzw_tree.iter->lzw_codes[*pixel_iter] = encoder->next_available_lzw_code;
            encoder->next_available_lzw_code += 1;

            // When the table is full, the encoder can choose to use the table as is, making no
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

                // Empty the LZW tree.

                LzwTreeNode *root = encoder->lzw_tree.root;
                memset(root->children, 0, sizeof(root->children));

                if (root != encoder->lzw_tree.last) {
                    encoder->lzw_tree.last->next = encoder->free_lzw_nodes;
                    encoder->lzw_tree.last = root;

                    encoder->free_lzw_nodes = root->next;
                    root->next = NULL;
                }

                encoder->lzw_code_bit_length = encoder->starting_lzw_code_bit_length;
                encoder->next_available_lzw_code = encoder->starting_lzw_code;
            }

            // Restart matching color sequences from the root of the tree.
            encoder->lzw_tree.iter = encoder->lzw_tree.root;
            encoder->next_lzw_code = LZW_CODE_INVALID;
        }
    }

    return pixel_iter - pixels;
}

void gif_encoder_finish_frame(GifEncoder *encoder, GifOutputBuffer *out_buffer) {
    assert(encoder->state == GIF_ENCODER_FRAME_STARTED);
    assert(gif_out_buffer_capacity_left(out_buffer) >= GIF_OUT_BUFFER_MIN_CAPACITY);

    // Encode the leftover sequence.
    if (encoder->lzw_tree.iter != encoder->lzw_tree.root) {
        gif_encoder_write_lzw_code(
            encoder,
            encoder->next_lzw_code,
            encoder->lzw_code_bit_length,
            out_buffer
        );
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

        isize pixels_consumed = gif_encoder_feed_frame(
            encoder,
            pixel_iter,
            pixels_end - pixel_iter,
            out_buffer
        );
        if (pixels_consumed >= 0) {
            pixel_iter += pixels_consumed;
        } else {
            return false;
        }
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
