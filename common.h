#ifndef COMMON_H
#define COMMON_H

#include <stdint.h> // uint8_t, uint16_t, uint32_t, int32_t
#include <stddef.h> // ptrdiff_t

#define UNUSED(x) ((void) (x))

typedef uint8_t u8;
typedef uint16_t u16;
typedef int32_t i32;
typedef uint32_t u32;
typedef ptrdiff_t isize;
typedef float f32;

#define FMT_ISIZE "%td"

#define sizeof(expr) ((isize) sizeof(expr))
#define countof(expr) (sizeof(expr) / sizeof((expr)[0]))

#endif // COMMON_H
