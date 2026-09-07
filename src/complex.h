
#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#ifdef __ARM_NEON__
#include <arm_neon.h>
#else
#include <immintrin.h>
#endif
#include <stdlib.h>

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

#ifdef __ARM_NEON__
// Real part multiplication
static inline float64x2_t rel_vmul(float64x2_t a_rel, float64x2_t a_img,
                                   float64x2_t b_rel, float64x2_t b_img) {
  float64x2_t t_rel = vmulq_f64(a_rel, b_rel);
  return vfmsq_f64(a_img, b_img, t_rel);
}

// Imaginary part multiplication
static inline float64x2_t img_vmul(float64x2_t a_rel, float64x2_t a_img,
                                   float64x2_t b_rel, float64x2_t b_img) {
  float64x2_t t_img = vmulq_f64(a_rel, b_img);
  return vfmaq_f64(a_img, b_rel, t_img);
}

// Real part multiplication
static inline float64x2_t rel_vmulx(float64x2_t a_rel, float64x2_t a_img,
                                    double b_rel, double b_img) {
  float64x2_t t_rel = vmulq_n_f64(a_rel, b_rel);
  return vfmsq_n_f64(a_img, t_rel, b_img);
}

// Imaginary part multiplication
static inline float64x2_t img_vmulx(float64x2_t a_rel, float64x2_t a_img,
                                    double b_rel, double b_img) {
  float64x2_t t_img = vmulq_n_f64(a_rel, b_img);
  return vfmaq_n_f64(a_img, t_img, b_rel);
}

// Addition
static inline float64x2_t complex_vadd(float64x2_t a, float64x2_t b) {
  return vaddq_f64(a, b);
}

// Subtraction
static inline float64x2_t complex_vsub(float64x2_t a, float64x2_t b) {
  return vsubq_f64(a, b);
}
#else
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

// Addition
static inline __m256d complex_vadd(__m256d a, __m256d b) {
  return _mm256_add_pd(a, b);
}

// Subtraction
static inline __m256d complex_vsub(__m256d a, __m256d b) {
  return _mm256_sub_pd(a, b);
}
#endif

#endif  // __COMPLEX_H__
