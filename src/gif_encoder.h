// Conventions and assumptions:
//  - Color arrays have either 3 or 4 components per pixel: red, green, blue, alpha (in this order).
//  - Functions which create or allocate something will return NULL if you run out of memory.
//  - Use encoder functions only when there is enough free memory left in the output buffer
//  (see GIF_OUT_BUFFER_MIN_CAPACITY).
//
// Customizing allocation:
//  You can provide a custom hook which is called right before allocating memory from the arena.
//  Provide the function name at compile time with the GIF_LIB_REPORT_ALLOC #define.
//  Hook signature: void report_alloc(void *arena, isize size)
//
// An output buffer exists outside of the encoder. This way the user can decide, if they want to use
// a dynamically growing buffer and gradually get the whole GIF into memory, or if they want to use
// a fixed buffer and flush it into the file every time it gets full.
//
// https://nullprogram.com/blog/2023/09/27
// I liked this guy's idea of using arenas in the library interface as opposed to abstracting
// allocation behind an ad hoc VMT malloc/realloc/free interface, and decided to try it out myself.
//
// An arena itself is just a pair of pointers: struct { u8 *begin; u8 *end; }
// No arena typedef is provided in this header, so that you could define it yourself and use it not
// only for GIF encoding, but anywhere else in your code.

#ifndef GIF_ENCODER_H
#define GIF_ENCODER_H

#include <stdbool.h>

// Redefinition of typedefs is a C11 feature.
// This is the official™ guard, which is used across different headers to protect u8 and friends.
// (Or just add a #define before including this header, if you already have short names defined.)
#ifndef SHORT_NAMES_FOR_PRIMITIVE_TYPES_WERE_DEFINED
    #define SHORT_NAMES_FOR_PRIMITIVE_TYPES_WERE_DEFINED
    #include <stdint.h>
    #include <stddef.h>

    typedef uint8_t   u8; typedef int8_t   i8;
    typedef uint16_t u16; typedef int16_t i16;
    typedef uint32_t u32; typedef int32_t i32;
    typedef uint64_t u64; typedef int64_t i64;

    typedef size_t   usize; typedef ptrdiff_t isize;
    typedef uintptr_t uptr;

    typedef float f32; typedef double f64;
#endif


// Every time you use encoder you need to have at least this amount of free memory left in buffer.
#define GIF_OUT_BUFFER_MIN_CAPACITY 1024

typedef struct {
    u8 *data;
    isize capacity;

    isize encoded_size;
    isize byte_pos;
    int bit_pos;
} GifOutputBuffer;

// Check if the buffer is full, read "encoded_size" bytes of encoded data from the "data" array...
isize gif_out_buffer_capacity_left(GifOutputBuffer const *out_buffer);
// ...and then reset the buffer to release space for the newly encoded data.
void gif_out_buffer_reset(GifOutputBuffer *out_buffer);

// These are just for convenience, you can manage the output buffer memory with malloc/realloc/free.
GifOutputBuffer *gif_out_buffer_create(isize min_capacity, void *arena);
bool gif_out_buffer_grow(GifOutputBuffer *out_buffer, isize min_capacity, void *arena);


#define GIF_MAX_WIDTH 0xffff
#define GIF_MAX_HEIGHT 0xffff
#define GIF_MAX_COLORS 256

typedef struct GifEncoder GifEncoder;

GifEncoder *gif_encoder_create(void *arena);

void gif_encoder_start(
    GifEncoder *encoder,
    isize width, isize height,
    u8 const *global_colors, isize global_color_count, int components,
    GifOutputBuffer *out_buffer
);

void gif_encoder_finish(GifEncoder *encoder, GifOutputBuffer *out_buffer);


// Call this before starting the frame. The default frame delay is 1 second.
void gif_frame_delay(GifEncoder *encoder, f32 seconds);

#define GIF_COLOR_INDEX_MAX 255
typedef u8 GifColorIndex;

void gif_encoder_start_frame(
    GifEncoder *encoder,
    u8 const *local_colors, isize local_color_count, int components,
    GifOutputBuffer *out_buffer
);

// Returns the number of pixels consumed. Call it repeatedly until all frame pixels are consumed.
// Returns -1 if encoder runs out of memory. You can't recover from this.
isize gif_encoder_feed_frame(
    GifEncoder *encoder,
    GifColorIndex const *pixels, isize pixel_count,
    GifOutputBuffer *out_buffer
);

void gif_encoder_finish_frame(GifEncoder *encoder, GifOutputBuffer *out_buffer);

// A shorthand for: start_frame -> feed_frame ... feed_frame -> finish_frame
//
// Fails and returns false if there is not enough memory in the buffer for the entire encoded frame
// or if we run out of memory in the process.
bool gif_encode_whole_frame(
    GifEncoder *encoder,
    u8 const *local_colors, isize color_count, int components,
    GifColorIndex *pixels,
    GifOutputBuffer *out_buffer
);


// Some predefined fixed color sRGB palettes.

u8 *srgb_palette_black_and_white(isize *color_count, int components, void *arena);
u8 *srgb_palette_monochrome(isize *color_count, int components, void *arena);
u8 *srgb_palette_web_safe(isize *color_count, int components, void *arena);


// Use these to convert the input image or sRGB palettes into the desired color space.
// Conversions back into sRGB are needed to put generated color tables into the GIF.

f32 *srgb_to_float(u8  const *srgb_colors, isize color_count, int components, void *arena);
u8  *float_to_srgb(f32 const *srgb_colors, isize color_count, int components, void *arena);

f32 *srgb_to_linear(u8  const *srgb_colors,   isize color_count, int components, void *arena);
u8  *linear_to_srgb(f32 const *linear_colors, isize color_count, int components, void *arena);

f32 *srgb_to_lab(u8  const *srgb_colors, isize color_count, int components, void *arena);
u8  *lab_to_srgb(f32 const *lab_colors,  isize color_count, int components, void *arena);

f32 *srgb_to_oklab(u8  const *srgb_colors,  isize color_count, int components, void *arena);
u8  *oklab_to_srgb(f32 const *oklab_colors, isize color_count, int components, void *arena);


// Generate a custom color palette that best fits the given image.
// Less colors passed to the palette generation algorithm the faster it works. Palette quality will
// differ depending on the specific algorithm.

u8 *colors_unique(
    u8 const *colors, isize color_count, int components,
    isize *unique_color_count,
    void *arena
);

// Writes the result directly into the "colors" array. Still needs extra memory for the hash set.
u8 *colors_unique_inplace(
    u8 *colors, isize color_count, int components,
    isize *unique_color_count,
    void *arena
);

f32 *palette_by_median_cut(
    f32 const *colors, isize color_count, int components,
    isize target_color_count,
    isize *colors_generated,
    void *arena
);

f32 *palette_by_k_means(
    f32 const *colors, isize color_count, int components,
    isize target_color_count,
    isize *colors_generated,
    bool k_means_plus_plus,
    void *arena
);

// Modified median-cut algorithm implemented following the description from here:
//  - http://leptonica.org/color-quantization.html
//  - http://leptonica.org/papers/mediancut.pdf
//  - https://tpgit.github.io/Leptonica/colorquant2_8c.html
u8 *palette_by_modified_median_cut(
    u8 const *colors, isize color_count, int components,
    isize target_color_count,
    isize *colors_generated,
    void *arena
);

// Octree palette generation implemented following these articles:
//  - https://www.cubic.org/docs/octree.htm
//  - https://observablehq.com/@tmcw/octree-color-quantization
u8 *palette_by_octree(
    u8 const *colors, isize color_count, int components,
    isize target_color_count,
    isize *colors_generated,
    void *arena
);


// Prepare an image to be fed into the GIF encoder.

GifColorIndex *image_quantize(
    f32 const *pixels, isize pixel_count, int components,
    f32 const *colors, isize color_count,
    void *arena
);

GifColorIndex *image_floyd_steinberg_dither(
    f32 const *image, isize width, isize height, int components,
    f32 const *colors, isize color_count,
    void *arena
);

GifColorIndex *image_ordered_dither(
    f32 const *image, isize width, isize height, int components,
    f32 const *colors, isize color_count,
    void *arena
);

#endif // GIF_ENCODER_H
