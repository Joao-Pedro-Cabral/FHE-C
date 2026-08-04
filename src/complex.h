
#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <stdint.h>

#include "fixed.h"

// Real part multiplication
static inline int64_t rel_mul(int64_t a_rel, int64_t a_img, int64_t b_rel,
                              int64_t b_img) {
  return fixed_sub(fixed_mul(a_rel, b_rel), fixed_mul(a_img, b_img));
}

// Imaginary part multiplication
static inline int64_t img_mul(int64_t a_rel, int64_t a_img, int64_t b_rel,
                              int64_t b_img) {
  return fixed_add(fixed_mul(a_rel, b_img), fixed_mul(a_img, b_rel));
}

#endif  // __COMPLEX_H__
