
#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <immintrin.h>

// We do static inline here for remove function call overhead

// Real part multiplication
static inline double rel_mul(double a_rel, double a_img, double b_rel,
                             double b_img) {
  return a_rel * b_rel - a_img * b_img;
}

// Imaginary part multiplication
static inline double img_mul(double a_rel, double a_img, double b_rel,
                             double b_img) {
  return a_rel * b_img + a_img * b_rel;
}

// Real part multiplication
static inline __m256d rel_vmul(__m256d a_rel, __m256d a_img, __m256d b_rel,
                               __m256d b_img) {
  __m256d t_rel = _mm256_mul_pd(a_img, b_img);
  return _mm256_fmsub_pd(a_rel, b_rel, t_rel);
}

// Imaginary part multiplication
static inline __m256d img_vmul(__m256d a_rel, __m256d a_img, __m256d b_rel,
                               __m256d b_img) {
  __m256d t_img = _mm256_mul_pd(a_rel, b_img);
  return _mm256_fmadd_pd(a_img, b_rel, t_img);
}

#endif  // __COMPLEX_H__
