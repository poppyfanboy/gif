#ifndef GIF_ENCODER_H
#define GIF_ENCODER_H

#include <stdbool.h>

#define GIF_MAX_WIDTH 0xffff
#define GIF_MAX_HEIGHT 0xffff
#define GIF_MAX_COLORS 256

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

// I saw this blog post:
// https://nullprogram.com/blog/2023/09/27
// and decided to try out the idea of explicitly using arenas in the library interface, as opposed
// to abstracting allocation behind a configurable malloc/realloc/free interface.
//
// I think that arena allocations should be made as explicit as possible, and because of that:
//
//  - Output buffer has to exist outside of the encoder; library user is responsible for resizing
//    the buffer, if needed. This way the user can decide, if they want to use a dynamic buffer and
//    gradually get the whole GIF into the memory, or if they want to use a fixed buffer and flush
//    it into the file from time to time.
//
//  - Encoder should not hold onto the pointer to the arena, it needs to preallocate whatever memory
//    it will need for the internal bookkeeping in advance.
//
// An arena itself is just a pair of pointers. I don't provide the typedef in this header, so that
// you had an option to define it yourself in your code.
//
//  typedef struct {
//      u8 *begin;
//      u8 *end;
//  } Arena;

// Color arrays here are assumed to have 3 components per pixel: red, green, blue (in this order).

// Some predefined fixed color palettes.
u8 *srgb_palette_black_and_white(isize *color_count, void *arena);
u8 *srgb_palette_monochrome(isize *color_count, void *arena);
u8 *srgb_palette_web_safe(isize *color_count, void *arena);

// Use these to convert the input image or sRGB palettes. Conversions back into sRGB might be needed
// to encode global/local color tables into the GIF.

f32 *srgb_to_float(u8 const *srgb_colors, isize color_count, void *arena);
u8 *float_to_srgb(f32 const *srgb_colors, isize color_count, void *arena);

f32 *srgb_to_linear(u8 const *srgb_colors, isize color_count, void *arena);
u8 *linear_to_srgb(f32 const *linear_colors, isize color_count, void *arena);

f32 *srgb_to_lab(u8 const *srgb_colors, isize color_count, void *arena);
u8 *lab_to_srgb(f32 const *lab_colors, isize color_count, void *arena);

f32 *srgb_to_oklab(u8 const *srgb_colors, isize color_count, void *arena);
u8 *oklab_to_srgb(f32 const *oklab_colors, isize color_count, void *arena);

// Generate a custom color palette that best fits the given image.
f32 *palette_by_median_cut(
    f32 const *pixels, isize pixel_count,
    isize target_color_count,
    isize *color_count,
    void *arena
);

// Prepare an image to be fed into the GIF encoder.

#define GIF_COLOR_INDEX_MAX 255
typedef u8 GifColorIndex;

GifColorIndex *image_quantize_for_gif(
    f32 const *pixels, isize pixel_count,
    f32 const *colors, isize color_count,
    void *arena
);

// The absolute minimum amount of memory needed to encode the whole GIF, given that you flush and
// reset the buffer after each call to any encoder function that writes something to the buffer.
#define GIF_OUT_BUFFER_MIN_CAPACITY 1024

typedef struct {
    u8 *data;
    isize encoded_size;
    isize byte_pos;
    int bit_pos;
    isize capacity;
} GifOutputBuffer;

// These make sure that any encoder function that writes something into the buffer can be called
// successfully at least once.
GifOutputBuffer *gif_out_buffer_create(isize min_capacity, void *arena);
isize gif_out_buffer_capacity_left(GifOutputBuffer const *out_buffer);
bool gif_out_buffer_grow(GifOutputBuffer *out_buffer, isize min_capacity, void *arena);

// Reset the buffer to free memory for new encoded data.
void gif_out_buffer_reset(GifOutputBuffer *out_buffer);

typedef struct GifEncoder GifEncoder;

isize gif_encoder_required_memory(void);
GifEncoder *gif_encoder_create(void *arena);

// Almost every encoder function below causes a buffer overflow, if there is not enough space in the
// buffer for it to successfully finish writing their portion of encoded data. You are supposed to
// use gif_out_buffer_capacity_left to check if there is enough memory (compare it against
// GIF_OUT_BUFFER_MIN_CAPACITY). I made it this way to reduce the number of possible internal states
// of the encoder.

void gif_encoder_start(
    GifEncoder *encoder,
    isize width,
    isize height,
    u8 const *global_colors,
    isize global_color_count,
    GifOutputBuffer *out_buffer
);

void gif_frame_delay(GifEncoder *encoder, f32 seconds);

void gif_encoder_start_frame(
    GifEncoder *encoder,
    u8 const *local_colors,
    isize local_color_count,
    GifOutputBuffer *out_buffer
);

// Reads a portion of the input image and returns the number of consumed pixels.
// Writes encoded data into the buffer.
// Call it repeatedly, until the whole frame is encoded.
isize gif_encoder_feed_frame(
    GifEncoder *encoder,
    GifColorIndex const *pixels,
    isize pixel_count,
    GifOutputBuffer *out_buffer
);

void gif_encoder_finish_frame(GifEncoder *encoder, GifOutputBuffer *out_buffer);

// Encodes the whole frame in a single function call, but gracefully fails and doesn't encode
// anything at all, if there isn't enough memory in the buffer. You are not supposed to continue
// execution as normal in case of failure: should have allocated more memory from the beginning or
// bought more RAM.
bool gif_encode_whole_frame(
    GifEncoder *encoder,
    u8 const *local_colors,
    isize color_count,
    GifColorIndex *pixels,
    GifOutputBuffer *out_buffer
);

void gif_encoder_finish(GifEncoder *encoder, GifOutputBuffer *out_buffer);

#endif // GIF_ENCODER_H
