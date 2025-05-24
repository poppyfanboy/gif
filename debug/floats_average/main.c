// https://orlp.net/blog/taming-float-sums
//
// More approaches:
// https://nhigham.com/2020/07/07/what-is-stochastic-rounding
// https://gitlab.com/radfordneal/xsum
//
// I'm probably doing something (or alotathing) wrong here:
// For some reason the average produced using Kahan summation is very often the exact same produced
// by the fsum ported from Python which looks extremely suspicious. And also average_pairwise
// usually produces less accurate results than average_block_pairwise, even though the latter is
// supposed to be a more performant but less accurate variation of the algorithm.
//
// But anyway, the conclusion is that for summing up f32s I just need to use f64 and it's going to
// be more than good enough for averaging colors. Splitting the single sum accumulator into multiple
// ones might improve things (compilers can't do this themselves without -fassociative-math or
// -ffast-math or whatever), but not by a lot.

#include <assert.h> // assert
#include <stdlib.h> // abort, malloc, qsort
#include <stddef.h> // size_t, NULL
#include <stdio.h>  // printf
#include <string.h> // memcpy
#include <math.h>   // fabsf, fabs, copysignf

#include <time.h>   // struct timespec, clock_gettime

#include <immintrin.h>

#if defined(__GNUC__) || defined(__clang__)

// https://github.com/google/benchmark/blob/main/docs/user_guide.md#preventing-optimization
// https://theunixzoo.co.uk/blog/2021-10-14-preventing-optimisations.html
#define DO_NOT_OPTIMIZE(value) asm volatile("" : "+m,r"(value) : : "memory")
#define CLOBBER_MEMORY() asm volatile("" : : : "memory")

#elif defined(_MSC_VER)

// https://github.com/google/benchmark/blob/176ad6c20c2a3a4a2613e041b90b76b4102aae3d/include/benchmark/benchmark.h#L614

#include <intrin.h> // _ReadWriteBarrier

volatile void const *volatile msvc_do_not_optimize_side_effect;

void msvc_do_not_optimize(char const *value) {
    msvc_do_not_optimize_side_effect = (void const *volatile)value;
    _ReadWriteBarrier();
}

#define DO_NOT_OPTIMIZE(value) msvc_do_not_optimize((char const *)&(value))
#define CLOBBER_MEMORY() _ReadWriteBarrier()

#else

#define DO_NOT_OPTIMIZE(value)
#define CLOBBER_MEMORY()

#endif

// Redefinition of typedefs is a C11 feature.
// This is the official™ guard, which is used across different headers to protect u8 and friends.
// (Or just add a #define before including this header, if you already have short names defined.)
#ifndef SHORT_NAMES_FOR_PRIMITIVE_TYPES_WERE_DEFINED
    #define SHORT_NAMES_FOR_PRIMITIVE_TYPES_WERE_DEFINED

    #include <stdint.h>
    #include <stddef.h>

    typedef int8_t i8;
    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef int16_t i16;
    typedef uint32_t u32;
    typedef int32_t i32;
    typedef uint64_t u64;
    typedef int64_t i64;

    typedef uintptr_t uptr;
    typedef size_t usize;
    typedef ptrdiff_t isize;

    typedef float f32;
    typedef double f64;
#endif

#define sizeof(expression) (isize)sizeof(expression)
#define lengthof(string) (sizeof(string) - 1)
#define countof(expression) (sizeof(expression) / sizeof((expression)[0]))

#define UNUSED(x) (void)x

int f32_compare(void const *left_ptr, void const *right_ptr) {
    f32 left = *(f32 *)left_ptr;
    f32 right = *(f32 *)right_ptr;

    return (int)copysignf(1.0F, left - right);
}

int f32_compare_decreasing(void const *left_ptr, void const *right_ptr) {
    f32 left = *(f32 *)left_ptr;
    f32 right = *(f32 *)right_ptr;

    return (int)copysignf(1.0F, right - left);
}

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct {
    u64 state;
    u64 inc;
} PCG32;

u32 pcg32_random(PCG32 *rng) {
    u64 oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    u32 xorshifted = (u32)(((oldstate >> 18U) ^ oldstate) >> 27U);
    u32 rot = oldstate >> 59U;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void pcg32_init(PCG32 *rng, u64 seed) {
    rng->state = 0U;
    rng->inc = 1U;
    pcg32_random(rng);
    rng->state += seed;
    pcg32_random(rng);
}

typedef struct {
    u8 *begin;
    u8 *end;
} Arena;

void *arena_alloc(Arena *arena, isize size) {
    assert(size > 0);

    // Padding for 16-byte memory address alignment.
    isize padding = (~(uptr)arena->begin + 1) & 0xff;
    isize memory_left = arena->end - arena->begin - padding;
    if (memory_left < 0 || memory_left < size) {
        abort();
    }

    void *ptr = arena->begin + padding;
    arena->begin += padding + size;
    return ptr;
}

f32 average_exact(f32 const *numbers, isize number_count, Arena arena);

typedef f32 (*BenchmarkProcedure)(f32 const *numbers, isize number_count, Arena arena);

void benchmark(BenchmarkProcedure procedure, u64 seed, Arena arena) {
    isize number_count = 100000000;
    f32 numbers_from = -1e9F * 0.5;
    f32 numbers_to = 1e9F * 1.5;

    PCG32 rng;
    pcg32_init(&rng, seed);

    isize accuracy_test_count = 10;
    u64 *accuracy_test_seeds = arena_alloc(&arena, accuracy_test_count * sizeof(u64));
    for (isize i = 0; i < accuracy_test_count; i += 1) {
        accuracy_test_seeds[i] = pcg32_random(&rng);
    }

    // Accuracy tests

    f64 accumulated_absolute_error = 0.0;
    for (isize accuracy_test = 0; accuracy_test < accuracy_test_count; accuracy_test += 1) {
        Arena arena_reset = arena;

        PCG32 rng_local;
        pcg32_init(&rng_local, accuracy_test_seeds[accuracy_test]);

        f32 *numbers = arena_alloc(&arena, number_count * sizeof(f32));
        for (isize i = 0; i < number_count; i += 1) {
            f32 random = (f32)(pcg32_random(&rng) >> 8) * 0x1.0p-24F;
            numbers[i] = random * (numbers_to - numbers_from) + numbers_from;
        }

        // Assume that the RNG is uniform and calcluate the exact value based on that instead?
        f32 exact_value = average_exact(numbers, number_count, arena);
        f32 prediction = procedure(numbers, number_count, arena);
        accumulated_absolute_error += fabs((f64)exact_value - (f64)prediction);

        arena = arena_reset;
    }

    f64 mean_absolute_error = accumulated_absolute_error / (f64)accuracy_test_count;
    printf("Mean absolute error: %.17g\n", mean_absolute_error);

    // Performance tests
    // I stole this from here (FormatJSONBenchmark.pas):
    // https://gitlab.com/freepascal.org/fpc/source/-/merge_requests/1018

    f32 *numbers = arena_alloc(&arena, number_count * sizeof(f32));
    for (isize i = 0; i < number_count; i += 1) {
        f32 random = (f32)(pcg32_random(&rng) >> 8) * 0x1.0p-24F;
            numbers[i] = random * (numbers_to - numbers_from) + numbers_from;
    }

    isize amplify = 10;

    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);

    isize repetitions = 0;
    f64 time;
    do {
        repetitions += 1;
        for (isize i = 0; i < amplify; i += 1) {
            f32 result = procedure(numbers, number_count, arena);
            DO_NOT_OPTIMIZE(result);
        }

        // Do I really need nanosecond-ish precision here?
        // I've heard that gettimeofday is a little bit faster, because it uses a different timer?
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        time = (f64)(now.tv_sec - start.tv_sec) + (f64)(now.tv_nsec - start.tv_nsec) * 1e-9;
    } while (time < 0.5);

    time = time / (f64)repetitions / (f64)amplify;

    char const *units[] = {"secs", "millis", "micros", "nanos"};
    isize unit_index = 0;
    while (unit_index < countof(units) - 1 && time < 0.9) {
        unit_index += 1;
        time *= 1e3;
    }

    printf("%.1f %s per call\n", time, units[unit_index]);
}

#define BENCHMARK(function, seed, arena)    \
    do {                                    \
        printf("%s\n", #function);          \
        benchmark(function, seed, arena);   \
        printf("\n");                       \
    } while (0)

f32 average_naive(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    f32 sum = 0.0F;

    for (isize i = 0; i < number_count; i += 1) {
        // At first I thought that doing
        // sum += numbers[i] / (f32)number_count
        // instead might be a better idea, but:
        //  - It takes more time, because now you need to do additional number_count divisions.
        //  - It might be even less accurate, because you're adding an even smaller number
        //  (numbers[i] / number_count) to the large number (sum).
        sum += numbers[i];
    }

    return sum / (f32)number_count;
}

// This is way too slow and probably does not do much to combat the error.
f32 average_naive_sorted(f32 const *numbers, isize number_count, Arena arena) {
    f32 *numbers_sorted = arena_alloc(&arena, number_count * sizeof(f32));
    memcpy(numbers_sorted, numbers, (size_t)(number_count * sizeof(f32)));
    qsort(numbers_sorted, (size_t)number_count, sizeof(f32), f32_compare);

    return average_naive(numbers_sorted, number_count, arena);
}

// I think this is it, just use f64 for computing the sum. Results are almost exact for reasonable
// input values. This whole thing was just a waste of time, I already do this anyway.
f32 average_naive_in_f64(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    f64 sum = 0.0;

    for (isize i = 0; i < number_count; i += 1) {
        sum += (f64)numbers[i];
    }

    return (f32)(sum / (f64)number_count);
}

f32 average_multi_x4_in_f64(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    f64 sum[4] = {0.0, 0.0, 0.0, 0.0};

    isize i = 0;
    for (i = 0; i < number_count; i += 4) {
        sum[0] += (f64)numbers[i];
        sum[1] += (f64)numbers[i + 1];
        sum[2] += (f64)numbers[i + 2];
        sum[3] += (f64)numbers[i + 3];
    }

    sum[0] += i < number_count ? (f64)numbers[i++] : 0.0;
    sum[1] += i < number_count ? (f64)numbers[i++] : 0.0;
    sum[2] += i < number_count ? (f64)numbers[i++] : 0.0;

    return (f32)((sum[0] + sum[1] + sum[2] + sum[3]) / (f64)number_count);
}

// I just wanted to try out using SSE intrinsics for the first time.
f32 average_sse2_in_f64(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    __m128d sum[2] = {_mm_setzero_pd(), _mm_setzero_pd()};

    isize i;
    for (i = 0; i < number_count; i += 4) {
        sum[0] = _mm_add_pd(sum[0], _mm_set_pd((f64)numbers[i], (f64)numbers[i + 1]));
        sum[1] = _mm_add_pd(sum[1], _mm_set_pd((f64)numbers[i + 2], (f64)numbers[i + 3]));
    }

    f64 tail[3];
    tail[0] = i < number_count ? (f64)numbers[i++] : 0.0;
    tail[1] = i < number_count ? (f64)numbers[i++] : 0.0;
    tail[2] = i < number_count ? (f64)numbers[i++] : 0.0;

    sum[0] = _mm_add_pd(sum[0], _mm_set_pd(tail[0], tail[1]));
    sum[1] = _mm_add_sd(sum[1], _mm_set_pd1(tail[2]));

    sum[0] = _mm_add_pd(sum[0], sum[1]);
    __m128d high = _mm_castps_pd(_mm_movehl_ps(_mm_castpd_ps(sum[0]), _mm_castpd_ps(sum[0])));
    sum[0] = _mm_add_sd(sum[0], high);

    return (f32)(_mm_cvtsd_f64(sum[0]) / (f64)number_count);
}

f32 sum_pairwise(f32 const *numbers, isize number_count) {
    if (number_count == 0) {
        return 0.0F;
    }
    if (number_count == 1) {
        return numbers[0];
    }

    isize split = number_count / 2;
    return sum_pairwise(numbers, split) + sum_pairwise(numbers + split, number_count - split);
}

f32 average_pairwise(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    return sum_pairwise(numbers, number_count) / (f32)number_count;
}

f32 sum_block_pairwise(f32 const *numbers, isize number_count) {
    isize block_size = 32;

    if (number_count <= block_size) {
        f32 sum = 0.0F;
        for (isize i = 0; i < number_count; i += 1) {
            sum += numbers[i];
        }
        return sum;
    } else {
        // Split at the next multiple of block_size:
        isize split = (number_count / 2 + block_size - 1) & ~(block_size - 1);
        return
            sum_block_pairwise(numbers, split) +
            sum_block_pairwise(numbers + split, number_count - split);
    }
}

f32 average_block_pairwise(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    return sum_block_pairwise(numbers, number_count) / (f32)number_count;
}

f32 average_kahan(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    f32 sum = 0.0F;
    f32 error = 0.0F;

    for (isize i = 0; i < number_count; i += 1) {
        f32 y = numbers[i] - error;
        f32 t = sum + y;
        error = (t - sum) - y;
        sum = t;
    }

    return sum / (f32)number_count;
}

isize isize_min(isize left, isize right) {
    return left < right ? left : right;
}

f32 average_block_kahan(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    isize block_size = 32;

    f32 sum = 0.0F;
    f32 error = 0.0F;

    for (isize i = 0; i < number_count; i += block_size) {
        f32 block_sum = 0.0F;
        for (isize j = i; j < isize_min(i + block_size, number_count); j += 1) {
            block_sum += numbers[j];
        }
        f32 y = block_sum - error;
        f32 t = sum + y;
        error = (t - sum) - y;
        sum = t;
    }

    return sum / (f32)number_count;
}

// https://github.com/python/cpython/blob/de19694cfbcaa1c85c3a4b7184a24ff21b1c0919/Modules/mathmodule.c#L1321
// https://en.wikipedia.org/wiki/2Sum
f32 average_python_fsum(f32 const *numbers, isize number_count, Arena arena) {
    struct {
        f32 *data;
        isize count;
        isize capacity;
    } partials;
    partials.data = arena_alloc(&arena, 32 * sizeof(f32));
    partials.count = 0;
    partials.capacity = 32;

    f32 const *number_iter = numbers;
    while (number_iter < numbers + number_count) {
        f32 x = *number_iter;

        f32 *partial_write_iter = partials.data;
        f32 const *partial_read_iter = partials.data;
        while (partial_read_iter < partials.data + partials.count) {
            f32 y = *partial_read_iter;
            if (fabsf(x) < fabsf(y)) {
                f32 swap = x;
                x = y;
                y = swap;
            }

            f32 sum = x + y;

            f32 error = y - (sum - x);
            if (error != 0.0F) {
                *partial_write_iter = error;
                partial_write_iter += 1;
            }

            x = sum;

            partial_read_iter += 1;
        }

        if (x != 0.0F) {
            if (partial_write_iter == partials.data + partials.capacity) {
                // × 1.375 growth factor
                isize new_capacity =
                    partials.capacity + (partials.capacity >> 2) + (partials.capacity >> 3);

                // realloc
                arena.begin += (new_capacity - partials.capacity) * sizeof(f32);
                if (arena.begin > arena.end) {
                    abort();
                }

                partials.capacity = new_capacity;
            }

            *partial_write_iter = x;
            partial_write_iter += 1;
            partials.count = partial_write_iter - partials.data;
        }

        number_iter += 1;
    }

    f32 sum = 0.0F;
    f32 const *partial_iter = partials.data + partials.count;

    if (partial_iter > partials.data) {
        partial_iter -= 1;
        sum = *partial_iter;

        f32 error = 0.0F;

        while (partial_iter > partials.data) {
            f32 x = sum;
            partial_iter -= 1;
            f32 y = *partial_iter;

            sum = x + y;
            f32 yr = sum - x;
            error = y - yr;
            if (error != 0.0F) {
                break;
            }
        }

        if (partial_iter > partials.data && (
            (error < 0.0F && *(partial_iter - 1) < 0.0F) ||
            (error > 0.0F && *(partial_iter - 1) > 0.0F)
        )) {
            f32 y = error * 2.0F;
            f32 x = sum + y;

            f32 yr = x - sum;
            if (y == yr) {
                sum = x;
            }
        }
    }

    return sum / (f32)number_count;
}

f32 average_exact(f32 const *numbers, isize number_count, Arena arena) {
    return average_python_fsum(numbers, number_count, arena);
}

int main(void) {
    isize arena_capacity = (isize)512 * 1024 * 1024;
    u8 *arena_memory = malloc((size_t)arena_capacity);
    Arena arena = {arena_memory, arena_memory + arena_capacity};

    u64 seed = 123;
    BENCHMARK(average_naive, seed, arena);
    BENCHMARK(average_naive_in_f64, seed, arena);
    BENCHMARK(average_kahan, seed, arena);
    BENCHMARK(average_block_kahan, seed, arena);
    BENCHMARK(average_multi_x4_in_f64, seed, arena);
    BENCHMARK(average_sse2_in_f64, seed, arena);
    BENCHMARK(average_pairwise, seed, arena);
    BENCHMARK(average_block_pairwise, seed, arena);
    BENCHMARK(average_python_fsum, seed, arena);

    return 0;
}
