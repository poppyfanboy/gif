// https://orlp.net/blog/taming-float-sums
//
// More approaches:
// https://nhigham.com/2020/07/07/what-is-stochastic-rounding
// https://gitlab.com/radfordneal/xsum
//
// The conclusion is that for computing the sum of f32s I just need to use f64. This approach is
// going to be more than good enough for averaging colors in the median-cut algorithm. Splitting a
// single sum accumulator variable into multiple separate accumulators improves performance as the
// compiler will now be able to properly vectorize the loop.

#include <assert.h> // assert
#include <stdlib.h> // abort, malloc, qsort
#include <stddef.h> // size_t, NULL
#include <stdio.h>  // printf, FILE, fopen, fscanf, feof
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

typedef enum {
    BENCHMARK_INPUT_DATA,
    BENCHMARK_INPUT_RANDOM,
} BenchmarkInputKind;

typedef struct {
    BenchmarkInputKind kind;
    isize performance_trials;
    union {
        struct {
            u64 seed;
            isize number_count;
            f32 numbers_from;
            f32 numbers_to;
            isize accuracy_trials;
        } random;

        struct {
            f32 *numbers;
            isize number_count;
        } data;
    } as;
} BenchmarkInput;

typedef f32 (*BenchmarkProcedure)(f32 const *numbers, isize number_count, Arena arena);

void benchmark(BenchmarkProcedure procedure, BenchmarkInput const *input, Arena arena) {
    // Accuracy tests

    f64 accumulated_absolute_error = 0.0;

    isize accuracy_trials;
    switch (input->kind) {
    case BENCHMARK_INPUT_DATA: {
        accuracy_trials = 1;
    } break;
    case BENCHMARK_INPUT_RANDOM: {
        accuracy_trials = input->as.random.accuracy_trials;
    } break;
    }

    if (input->kind == BENCHMARK_INPUT_DATA) {
        isize number_count = input->as.data.number_count;
        f32 *numbers = input->as.data.numbers;

        f32 expected = average_exact(numbers, number_count, arena);
        f32 actual = procedure(numbers, number_count, arena);

        printf("================\n");
        printf("expected value: %.9g\n", expected);
        printf("actual value: %.9g\n", actual);

        accumulated_absolute_error += fabs((f64)expected - (f64)actual);
    } else if (input->kind == BENCHMARK_INPUT_RANDOM) {
        PCG32 rng;
        pcg32_init(&rng, input->as.random.seed);

        isize number_count = input->as.random.number_count;

        u64 *accuracy_test_seeds = arena_alloc(&arena, accuracy_trials * sizeof(u64));
        for (isize i = 0; i < accuracy_trials; i += 1) {
            accuracy_test_seeds[i] = pcg32_random(&rng);
        }

        for (isize accuracy_test = 0; accuracy_test < accuracy_trials; accuracy_test += 1) {
            Arena arena_reset = arena;

            PCG32 rng_local;
            pcg32_init(&rng_local, accuracy_test_seeds[accuracy_test]);

            f32 numbers_from = input->as.random.numbers_from;
            f32 numbers_to = input->as.random.numbers_to;

            f32 *numbers = arena_alloc(&arena, number_count * sizeof(f32));
            for (isize i = 0; i < number_count; i += 1) {
                f32 random = (f32)(pcg32_random(&rng_local) >> 8) * 0x1.0p-24F;
                numbers[i] = random * (numbers_to - numbers_from) + numbers_from;
            }

            // Assume that the RNG is uniform and calcluate the exact value
            // as "(numbers_from + numbers_to) / 2"?
            //
            // Thinking about it once again, no, it's not even about potentially RNG not being perfectly
            // uniform, it's about the fact that we have a specific discrete case at hand instead of
            // some theoretical unreachable "number_count -> ∞" case.
            //
            // BTW is this how it's deduced? I'm not sure, it's not obvious for me for some reason.
            //
            // ∫0..1(from + (to - from) * x)dx =
            // from * x + (to - from) * x^2 / 2 | 0..1 =
            // from + (to - from) / 2 =
            // (from + to) / 2

            f32 exact_value = average_exact(numbers, number_count, arena);
            f32 prediction = procedure(numbers, number_count, arena);
            accumulated_absolute_error += fabs((f64)exact_value - (f64)prediction);

            arena = arena_reset;
        }
    }

    f64 mean_absolute_error = accumulated_absolute_error / (f64)accuracy_trials;
    printf("================\n");
    printf("Mean absolute error: %.17g\n", mean_absolute_error);

    // Performance tests
    // I stole this from here (FormatJSONBenchmark.pas):
    // https://gitlab.com/freepascal.org/fpc/source/-/merge_requests/1018

    isize performance_trials = input->performance_trials;

    isize number_count;
    f32 *numbers;
    switch (input->kind) {
    case BENCHMARK_INPUT_DATA: {
        number_count = input->as.data.number_count;
        numbers = input->as.data.numbers;
    } break;
    case BENCHMARK_INPUT_RANDOM: {
        PCG32 rng;
        pcg32_init(&rng, input->as.random.seed);

        number_count = input->as.random.number_count;
        numbers = arena_alloc(&arena, number_count * sizeof(f32));

        f32 numbers_from = input->as.random.numbers_from;
        f32 numbers_to = input->as.random.numbers_to;

        for (isize i = 0; i < number_count; i += 1) {
            f32 random = (f32)(pcg32_random(&rng) >> 8) * 0x1.0p-24F;
                numbers[i] = random * (numbers_to - numbers_from) + numbers_from;
        }
    } break;
    }

    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);

    isize repetitions = 0;
    f64 time;
    do {
        repetitions += 1;
        for (isize i = 0; i < performance_trials; i += 1) {
            f32 result = procedure(numbers, number_count, arena);
            DO_NOT_OPTIMIZE(result);
        }

        // Do I really need nanosecond-ish precision here?
        // I've heard that gettimeofday is a little bit faster, because it uses a different timer?
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        time = (f64)(now.tv_sec - start.tv_sec) + (f64)(now.tv_nsec - start.tv_nsec) * 1e-9;
    } while (time < 0.5);

    time = time / (f64)repetitions / (f64)performance_trials;

    char const *units[] = {"secs", "millis", "micros", "nanos"};
    isize unit_index = 0;
    while (unit_index < countof(units) - 1 && time < 0.9) {
        unit_index += 1;
        time *= 1e3;
    }

    printf("================\n");
    printf("%.1f %s per call\n", time, units[unit_index]);
}

#define BENCHMARK(function, input, arena)   \
    do {                                    \
        printf("%s\n", #function);          \
        benchmark(function, input, arena);  \
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

    // Cast to f64 in case number_count is so large that it cannot be represented in f32 with
    // reasonable precision?
    return (f32)((f64)sum / (f64)number_count);
}

// This is way too slow and probably does not do much to combat the error.
f32 average_naive_sorted(f32 const *numbers, isize number_count, Arena arena) {
    f32 *numbers_sorted = arena_alloc(&arena, number_count * sizeof(f32));
    memcpy(numbers_sorted, numbers, (size_t)(number_count * sizeof(f32)));
    qsort(numbers_sorted, (size_t)number_count, sizeof(f32), f32_compare);

    return average_naive(numbers_sorted, number_count, arena);
}

// I think this is it, just use f64 for computing the sum. Results are almost exact for reasonable
// input values. This whole thing was just a waste of time, I already use this approach anyway.
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

    f32 const *number_iter = numbers;
    f32 const *numbers_end = numbers + number_count;

    while (numbers_end - number_iter >= 4) {
        sum[0] += (f64)number_iter[0];
        sum[1] += (f64)number_iter[1];
        sum[2] += (f64)number_iter[2];
        sum[3] += (f64)number_iter[3];

        number_iter += 4;
    }

    if (number_iter < numbers_end) {
        sum[0] += (f64)*number_iter;
        number_iter += 1;
    }
    if (number_iter < numbers_end) {
        sum[1] += (f64)*number_iter;
        number_iter += 1;
    }
    if (number_iter < numbers_end) {
        sum[2] += (f64)*number_iter;
        number_iter += 1;
    }

    return (f32)((sum[0] + sum[1] + sum[2] + sum[3]) / (f64)number_count);
}

#ifdef AVX

// Pretty much useless. Compilers already do a good job of vectorizing average_multi_x4_in_f64.
f32 average_multi_x16_avx_in_f64(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    __m256d sum[4] = {
        _mm256_setzero_pd(),
        _mm256_setzero_pd(),
        _mm256_setzero_pd(),
        _mm256_setzero_pd(),
    };

    f32 const *number_iter = numbers;
    f32 const *numbers_end = numbers + number_count;

    while (numbers_end - number_iter >= 16) {
        sum[0] = _mm256_add_pd(sum[0], _mm256_cvtps_pd(_mm_loadu_ps(number_iter +  0)));
        sum[1] = _mm256_add_pd(sum[1], _mm256_cvtps_pd(_mm_loadu_ps(number_iter +  4)));
        sum[2] = _mm256_add_pd(sum[2], _mm256_cvtps_pd(_mm_loadu_ps(number_iter +  8)));
        sum[3] = _mm256_add_pd(sum[3], _mm256_cvtps_pd(_mm_loadu_ps(number_iter + 12)));

        number_iter += 16;
    }

    f64 tail[16];
    for (isize i = 0; i < 16; i += 1) {
        if (number_iter < numbers_end) {
            tail[i] = (f64)*number_iter;
            number_iter += 1;
        } else {
            tail[i] = 0.0;
        }
    }
    sum[0] = _mm256_add_pd(sum[0], _mm256_loadu_pd(tail +  0));
    sum[1] = _mm256_add_pd(sum[1], _mm256_loadu_pd(tail +  4));
    sum[2] = _mm256_add_pd(sum[2], _mm256_loadu_pd(tail +  8));
    sum[3] = _mm256_add_pd(sum[3], _mm256_loadu_pd(tail + 12));

    sum[0] = _mm256_add_pd(sum[0], sum[1]);
    sum[2] = _mm256_add_pd(sum[2], sum[3]);
    sum[0] = _mm256_add_pd(sum[0], sum[2]);

    // hadd([a, b, c, d]) -> [a + b, a + b, c + d, c + d]
    sum[0] = _mm256_hadd_pd(sum[0], sum[0]);
    // high = [c + d, c + d]
    __m128d high = _mm256_extractf128_pd(sum[0], 1);
    // low = [a + b, a + b]
    __m128d low = _mm256_castpd256_pd128(sum[0]);

    __m128d result = _mm_add_sd(low, high);
    return (f32)(_mm_cvtsd_f64(result) / (f64)number_count);
}

#endif // AVX

f32 sum_pairwise(f32 const *numbers, isize number_count) {
    if (number_count == 0) {
        return 0.0F;
    }
    if (number_count == 1) {
        return numbers[0];
    }

    // Split at some power of two.
    //
    // Splitting in half often gives even worse precision than sum_block_pairwise, so try to split
    // at powers of two? This seems to give more precise results *sometimes*.
    isize isize_max = (isize)(~(usize)0 >> 1);
    isize split = (isize_max >> __builtin_clzll((usize)number_count - 1)) + 1;

    return sum_pairwise(numbers, split) + sum_pairwise(numbers + split, number_count - split);
}

f32 average_pairwise(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    return (f32)(sum_pairwise(numbers, number_count) / (f64)number_count);
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

    return (f32)(sum_block_pairwise(numbers, number_count) / (f64)number_count);
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

    return (f32)(sum / (f64)number_count);
}

isize isize_min(isize left, isize right) {
    return left < right ? left : right;
}

f32 average_block_kahan(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    isize const block_size = 32;

    f32 sum = 0.0F;
    f32 error = 0.0F;

    f32 const *number_iter = numbers;
    f32 const *numbers_end = numbers + number_count;

    while (numbers_end - number_iter >= block_size) {
        f32 block_sum = 0.0F;
        for (isize i = 0; i < block_size; i += 1) {
            block_sum += number_iter[i];
        }

        f32 y = block_sum - error;
        f32 t = sum + y;
        error = (t - sum) - y;
        sum = t;

        number_iter += block_size;
    }

    f32 tail_sum = 0.0F;
    while (number_iter < numbers_end) {
        tail_sum += *number_iter;
        number_iter += 1;
    }
    sum += tail_sum - error;

    return (f32)(sum / (f64)number_count);
}

#ifdef AVX

// This one computes mulitple Kahan sums simultaneously and combines them afterwards, so it's going
// to be less precise than average_block_kahan.
//
// I don't know if there is some clever way to combine multiple Kahan sums, so I'm just using Kahan
// once again to sum up the "sum" terms; as for the "error" terms: they just get discarded.
f32 average_block_kahan_avx(f32 const *numbers, isize number_count, Arena arena) {
    UNUSED(arena);

    __m256 sum = _mm256_setzero_ps();
    __m256 error = _mm256_setzero_ps();

    f32 const *number_iter = numbers;
    f32 const *numbers_end = numbers + number_count;

    while (numbers_end - number_iter >= 32) {
        __m256 block_sum[4];
        block_sum[0] = _mm256_loadu_ps(number_iter +  0);
        block_sum[1] = _mm256_loadu_ps(number_iter +  8);
        block_sum[2] = _mm256_loadu_ps(number_iter + 16);
        block_sum[3] = _mm256_loadu_ps(number_iter + 24);

        block_sum[0] = _mm256_add_ps(block_sum[0], block_sum[1]);
        block_sum[2] = _mm256_add_ps(block_sum[2], block_sum[3]);
        block_sum[0] = _mm256_add_ps(block_sum[0], block_sum[2]);

        __m256 y = _mm256_sub_ps(block_sum[0], error);
        __m256 t = _mm256_add_ps(sum, y);
        error = _mm256_sub_ps(_mm256_sub_ps(t, sum), y);
        sum = t;

        number_iter += 32;
    }
    while (numbers_end - number_iter >= 8) {
        __m256 y = _mm256_sub_ps(_mm256_loadu_ps(number_iter), error);
        __m256 t = _mm256_add_ps(sum, y);
        error = _mm256_sub_ps(_mm256_sub_ps(t, sum), y);
        sum = t;

        number_iter += 8;
    }

    f32 sum_scalar[8];
    _mm256_storeu_ps(sum_scalar, sum);

    // Handle the tail
    {
        f32 error_scalar[8];
        _mm256_storeu_ps(error_scalar, error);

        for (isize i = 0; i < numbers_end - number_iter; i += 1) {
            sum_scalar[i] += number_iter[i] - error_scalar[i];
        }
    }

    f32 final_sum = 0.0F;
    f32 final_error = 0.0F;
    for (isize j = 0; j < 8; j += 1) {
        f32 y = sum_scalar[j] - final_error;
        f32 t = final_sum + y;
        final_error = (t - final_sum) - y;
        final_sum = t;
    }

    return (f32)(final_sum / (f64)number_count);
}

#endif // AVX

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
        }
        partials.count = partial_write_iter - partials.data;

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

    return (f32)(sum / (f64)number_count);
}

f32 average_exact(f32 const *numbers, isize number_count, Arena arena) {
    return average_python_fsum(numbers, number_count, arena);
}

int main(void) {
    isize arena_capacity = (isize)512 * 1024 * 1024;
    u8 *arena_memory = malloc((size_t)arena_capacity);
    Arena arena = {arena_memory, arena_memory + arena_capacity};

#if 1
    BenchmarkInput input = {
        .kind = BENCHMARK_INPUT_RANDOM,
        .performance_trials = 10000,
        .as.random = {
            .seed = 123,

            // more numbers => larger sum => less precision difference between different methods
            .number_count = 1000000,
            .numbers_from = -1e6F * 0.5,
            .numbers_to = 1e6F * 1.5,

            .accuracy_trials = 200,
        },
    };
#elif 1
    BenchmarkInput input = {
        .kind = BENCHMARK_INPUT_DATA,
        .performance_trials = 1000000000,
        .as.data = {
            .numbers = (f32[]){1e25F, 3.0F, -1e25F},
            .number_count = 3,
        },
    };
#else
    FILE *file = fopen("input.txt", "rb");
    if (file == NULL) {
        return 1;
    }

    isize number_count = 0;
    isize capacity = 32;
    f32 *numbers = arena_alloc(&arena, capacity * sizeof(f32));

    while (feof(file) == 0) {
        if (number_count == capacity) {
            isize new_capacity = capacity * 2;
            arena.begin += (new_capacity - capacity) * sizeof(f32);
            if (arena.begin > arena.end) {
                abort();
            }
            capacity = new_capacity;
        }

        fscanf(file, "%f\n", &numbers[number_count]);
        number_count += 1;
    }

    BenchmarkInput input = {
        .kind = BENCHMARK_INPUT_DATA,
        .performance_trials = 1000000,
        .as.data = {
            .numbers = numbers,
            .number_count = number_count,
        },
    };
#endif

    BENCHMARK(average_naive, &input, arena);
    BENCHMARK(average_naive_in_f64, &input, arena);
    BENCHMARK(average_multi_x4_in_f64, &input, arena);
#ifdef AVX
    BENCHMARK(average_multi_x16_avx_in_f64, &input, arena);
#endif

    BENCHMARK(average_kahan, &input, arena);
    BENCHMARK(average_block_kahan, &input, arena);
#ifdef AVX
    BENCHMARK(average_block_kahan_avx, &input, arena);
#endif

    BENCHMARK(average_pairwise, &input, arena);
    BENCHMARK(average_block_pairwise, &input, arena);

    BENCHMARK(average_python_fsum, &input, arena);

    return 0;
}
