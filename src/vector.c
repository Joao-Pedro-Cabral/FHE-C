
#include "vector.h"

#ifdef __NEON_MODE__
#include <arm_neon.h>
#elif __AVX2_MODE__
#include <immintrin.h>
#endif
#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include "mod.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "complex.h"
#include "utils.h"

// m = n/2 (half cyclotomic size)
void vec_new(vec_t* v, uint32_t m) {
  v->rel = malloc(2 * m * sizeof(double));
  v->img = &v->rel[m];
}

void vec_free(vec_t* v) { free(v->rel); }

// m = n/2 (half cyclotomic size)
void vec_gen_roots(vec_t* w, uint32_t m) {
  // Position 0 is not used
  w->rel[0] = 0.0;
  w->img[0] = 0.0;
#pragma omp for
  for (uint32_t i = 1; i < m; i++) {
    uint32_t log2_i = floor_log2(i);
    uint32_t pow2_i = i - (1 << log2_i);
    uint32_t mod2_i = 1 << (log2_i + 3);
    uint32_t pow5 = mod_pow2(5, pow2_i, mod2_i);
    uint32_t gap = m >> (log2_i + 1);
    uint32_t idx = pow5 * gap;
    double cosw = cos((idx * M_PI) / (2 * m));
    double sinw = sin((idx * M_PI) / (2 * m));
    w->rel[i] = cosw;
    w->img[i] = sinw;
  }
}

#define CT_BUTTERFLY(vj_rel, vj_img, vk_rel, vk_img, vw_rel, vw_img, vt_rel, \
                     vt_img)                                                 \
  do {                                                                       \
    (vt_rel) = rel_vmul((vk_rel), (vk_img), (vw_rel), (vw_img));             \
    (vt_img) = img_vmul((vk_rel), (vk_img), (vw_rel), (vw_img));             \
    (vk_rel) = complex_vsub((vj_rel), (vt_rel));                             \
    (vk_img) = complex_vsub((vj_img), (vt_img));                             \
    (vj_rel) = complex_vadd((vj_rel), (vt_rel));                             \
    (vj_img) = complex_vadd((vj_img), (vt_img));                             \
  } while (0)

#ifdef __NEON_MODE__
#define CT_BUTTERFLY_X(vj_rel, vj_img, vk_rel, vk_img, w_rel, w_img, vt_rel, \
                       vt_img)                                               \
  do {                                                                       \
    (vt_rel) = rel_vmulx((vk_rel), (vk_img), (w_rel), (w_img));              \
    (vt_img) = img_vmulx((vk_rel), (vk_img), (w_rel), (w_img));              \
    (vk_rel) = complex_vsub((vj_rel), (vt_rel));                             \
    (vk_img) = complex_vsub((vj_img), (vt_img));                             \
    (vj_rel) = complex_vadd((vj_rel), (vt_rel));                             \
    (vj_img) = complex_vadd((vj_img), (vt_img));                             \
  } while (0)
#endif

// m = n/2 (half cyclotomic size)
void fft(vec_t* v, const vec_t* w, uint32_t m) {
#ifdef __NEON_MODE__
  if (m < 16) {
#pragma omp single
    for (uint32_t len = 1; len < m; len <<= 1) {
      for (uint32_t i = 0; i < m; i += 2 * len) {
        for (uint32_t j = i; j < i + len; j++) {
          double w_rel = w->rel[len + j - i];
          double w_img = w->img[len + j - i];
          double j_rel = v->rel[j];
          double j_img = v->img[j];
          double k_rel = v->rel[j + len];
          double k_img = v->img[j + len];
          double t_rel = rel_mul(k_rel, k_img, w_rel, w_img);
          double t_img = img_mul(k_rel, k_img, w_rel, w_img);
          v->rel[j + len] = j_rel - t_rel;
          v->img[j + len] = j_img - t_img;
          v->rel[j] = j_rel + t_rel;
          v->img[j] = j_img + t_img;
        }
      }
    }
    return;
  }

  // Based on the number of threads and the format of n,
  // we decide the almost-optimal block-division
  uint32_t n_threads = (uint32_t)omp_get_num_threads();
  uint32_t tid = (uint32_t)omp_get_thread_num();
  uint32_t chunk_size = m / n_threads;
  uint32_t log2_m = floor_log2(m / 16);
  uint32_t m3_rem = log2_m % 3;
  uint32_t m3 = 1 << (log2_m + 4 - m3_rem);
  bool sync_before_rem = chunk_size >= m3 && chunk_size < m;
  bool sync_last_3x = chunk_size < m3 && chunk_size < m;
  if (sync_before_rem)
    chunk_size = m3;
  else if (sync_last_3x)
    chunk_size = m3 >> 3;

  // Merge 4 layers using all the 32 ARM NEON (ARMv8) registers
  uint32_t k1 = 1;
  uint32_t k2 = 2;
  uint32_t k3 = 4;
  uint32_t k4 = 8;
  double w1_rel = w->rel[k1];
  double w1_img = w->img[k1];
  float64x2_t vw2_rel = vld1q_f64(&w->rel[k2]);
  float64x2_t vw2_img = vld1q_f64(&w->img[k2]);
  for (uint32_t jj = tid * chunk_size; jj < m; jj += n_threads * chunk_size) {
    for (uint32_t j = jj; j < jj + chunk_size; j += 16) {
      // Load
      // First Quarter
      float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
      float64x2_t vj0_img = vld1q_f64(&v->img[j]);
      float64x2_t vk0_rel = vld1q_f64(&v->rel[j + 2]);
      float64x2_t vk0_img = vld1q_f64(&v->img[j + 2]);
      // Second Quarter
      float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 4]);
      float64x2_t vj1_img = vld1q_f64(&v->img[j + 4]);
      float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 6]);
      float64x2_t vk1_img = vld1q_f64(&v->img[j + 6]);
      // Third Quarter
      float64x2_t vj2_rel = vld1q_f64(&v->rel[j + 8]);
      float64x2_t vj2_img = vld1q_f64(&v->img[j + 8]);
      float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 10]);
      float64x2_t vk2_img = vld1q_f64(&v->img[j + 10]);
      // Forth Quarter
      float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 12]);
      float64x2_t vj3_img = vld1q_f64(&v->img[j + 12]);
      float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 14]);
      float64x2_t vk3_img = vld1q_f64(&v->img[j + 14]);
      // Tranpose
      // First Quarter
      float64x2_t vt0_rel = vj0_rel;
      float64x2_t vt0_img = vj0_img;
      vj0_rel = vzip1q_f64(vt0_rel, vk0_rel);
      vj0_img = vzip1q_f64(vt0_img, vk0_img);
      vk0_rel = vzip2q_f64(vt0_rel, vk0_rel);
      vk0_img = vzip2q_f64(vt0_img, vk0_img);
      // Second Quarter
      float64x2_t vt1_rel = vj1_rel;
      float64x2_t vt1_img = vj1_img;
      vj1_rel = vzip1q_f64(vt1_rel, vk1_rel);
      vj1_img = vzip1q_f64(vt1_img, vk1_img);
      vk1_rel = vzip2q_f64(vt1_rel, vk1_rel);
      vk1_img = vzip2q_f64(vt1_img, vk1_img);
      // Third Quarter
      float64x2_t vt2_rel = vj2_rel;
      float64x2_t vt2_img = vj2_img;
      vj2_rel = vzip1q_f64(vt2_rel, vk2_rel);
      vj2_img = vzip1q_f64(vt2_img, vk2_img);
      vk2_rel = vzip2q_f64(vt2_rel, vk2_rel);
      vk2_img = vzip2q_f64(vt2_img, vk2_img);
      // Forth Quarter
      float64x2_t vt3_rel = vj3_rel;
      float64x2_t vt3_img = vj3_img;
      vj3_rel = vzip1q_f64(vt3_rel, vk3_rel);
      vj3_img = vzip1q_f64(vt3_img, vk3_img);
      vk3_rel = vzip2q_f64(vt3_rel, vk3_rel);
      vk3_img = vzip2q_f64(vt3_img, vk3_img);
      // First layer (len = 1), one root
      float64x2_t vw_rel;
      float64x2_t vw_img;
      // First Quarter
      CT_BUTTERFLY_X(vj0_rel, vj0_img, vk0_rel, vk0_img, w1_rel, w1_img,
                     vt0_rel, vt0_img);
      // Second Quarter
      CT_BUTTERFLY_X(vj1_rel, vj1_img, vk1_rel, vk1_img, w1_rel, w1_img,
                     vt1_rel, vt1_img);
      // Third Quarter
      CT_BUTTERFLY_X(vj2_rel, vj2_img, vk2_rel, vk2_img, w1_rel, w1_img,
                     vt2_rel, vt2_img);
      // Forth Quarter
      CT_BUTTERFLY_X(vj3_rel, vj3_img, vk3_rel, vk3_img, w1_rel, w1_img,
                     vt3_rel, vt3_img);
      // Transpose back
      // First Quarter
      vt0_rel = vj0_rel;
      vt0_img = vj0_img;
      vj0_rel = vzip1q_f64(vt0_rel, vk0_rel);
      vj0_img = vzip1q_f64(vt0_img, vk0_img);
      vk0_rel = vzip2q_f64(vt0_rel, vk0_rel);
      vk0_img = vzip2q_f64(vt0_img, vk0_img);
      // Second Quarter
      vt1_rel = vj1_rel;
      vt1_img = vj1_img;
      vj1_rel = vzip1q_f64(vt1_rel, vk1_rel);
      vj1_img = vzip1q_f64(vt1_img, vk1_img);
      vk1_rel = vzip2q_f64(vt1_rel, vk1_rel);
      vk1_img = vzip2q_f64(vt1_img, vk1_img);
      // Third Quarter
      vt2_rel = vj2_rel;
      vt2_img = vj2_img;
      vj2_rel = vzip1q_f64(vt2_rel, vk2_rel);
      vj2_img = vzip1q_f64(vt2_img, vk2_img);
      vk2_rel = vzip2q_f64(vt2_rel, vk2_rel);
      vk2_img = vzip2q_f64(vt2_img, vk2_img);
      // Forth Quarter
      vt3_rel = vj3_rel;
      vt3_img = vj3_img;
      vj3_rel = vzip1q_f64(vt3_rel, vk3_rel);
      vj3_img = vzip1q_f64(vt3_img, vk3_img);
      vk3_rel = vzip2q_f64(vt3_rel, vk3_rel);
      vk3_img = vzip2q_f64(vt3_img, vk3_img);
      // Second layer (len = 2), two roots
      // First Quarter
      CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw2_rel, vw2_img,
                   vt0_rel, vt0_img);
      // Second Quarter
      CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw2_rel, vw2_img,
                   vt1_rel, vt1_img);
      // Third Quarter
      CT_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw2_rel, vw2_img,
                   vt2_rel, vt2_img);
      // Forth Quarter
      CT_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw2_rel, vw2_img,
                   vt3_rel, vt3_img);
      // Third layer (len = 4), four roots
      // First Quarter
      vw_rel = vld1q_f64(&w->rel[k3]);
      vw_img = vld1q_f64(&w->img[k3]);
      CT_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw_rel, vw_img, vt0_rel,
                   vt0_img);
      // Third Quarter
      CT_BUTTERFLY(vj2_rel, vj2_img, vj3_rel, vj3_img, vw_rel, vw_img, vt2_rel,
                   vt2_img);
      // Second Quarter
      vw_rel = vld1q_f64(&w->rel[k3 + 2]);
      vw_img = vld1q_f64(&w->img[k3 + 2]);
      CT_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                   vt1_img);
      // Forth Quarter
      CT_BUTTERFLY(vk2_rel, vk2_img, vk3_rel, vk3_img, vw_rel, vw_img, vt3_rel,
                   vt3_img);
      // Forth layer (len = 8), eight roots
      // First Quarter
      vw_rel = vld1q_f64(&w->rel[k4]);
      vw_img = vld1q_f64(&w->img[k4]);
      CT_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw_rel, vw_img, vt0_rel,
                   vt0_img);
      // Second Quarter
      vw_rel = vld1q_f64(&w->rel[k4 + 2]);
      vw_img = vld1q_f64(&w->img[k4 + 2]);
      CT_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw_rel, vw_img, vt1_rel,
                   vt1_img);
      // Third Quarter
      vw_rel = vld1q_f64(&w->rel[k4 + 4]);
      vw_img = vld1q_f64(&w->img[k4 + 4]);
      CT_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw_rel, vw_img, vt2_rel,
                   vt2_img);
      // Quarter Quarter
      vw_rel = vld1q_f64(&w->rel[k4 + 6]);
      vw_img = vld1q_f64(&w->img[k4 + 6]);
      CT_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw_rel, vw_img, vt3_rel,
                   vt3_img);
      // Store
      // First Quarter
      vst1q_f64(&v->rel[j], vj0_rel);
      vst1q_f64(&v->img[j], vj0_img);
      vst1q_f64(&v->rel[j + 2], vk0_rel);
      vst1q_f64(&v->img[j + 2], vk0_img);
      // Second Quarter
      vst1q_f64(&v->rel[j + 4], vj1_rel);
      vst1q_f64(&v->img[j + 4], vj1_img);
      vst1q_f64(&v->rel[j + 6], vk1_rel);
      vst1q_f64(&v->img[j + 6], vk1_img);
      // Third Quarter
      vst1q_f64(&v->rel[j + 8], vj2_rel);
      vst1q_f64(&v->img[j + 8], vj2_img);
      vst1q_f64(&v->rel[j + 10], vk2_rel);
      vst1q_f64(&v->img[j + 10], vk2_img);
      // Forth Quarter
      vst1q_f64(&v->rel[j + 12], vj3_rel);
      vst1q_f64(&v->img[j + 12], vj3_img);
      vst1q_f64(&v->rel[j + 14], vk3_rel);
      vst1q_f64(&v->img[j + 14], vk3_img);
    }
  }

  // Merge 3 consecutive layers
  uint32_t jj_inc = n_threads * chunk_size;
  uint32_t klen_inc = 0;
  for (uint32_t len = 16; len < m3; len <<= 3) {
    uint32_t i_inc = 8 * len;
    if (sync_last_3x && len == (m3 >> 3)) {
      jj_inc = m3;
      chunk_size = chunk_size / n_threads;  // chunk_size = n3 >> 3
      i_inc = 8 * chunk_size;
      klen_inc = tid * chunk_size;
#pragma omp barrier
    }
    for (uint32_t jj = tid * chunk_size; jj < m; jj += jj_inc) {
      for (uint32_t i = jj; i < jj + chunk_size; i += i_inc) {
        uint32_t klen = len - i + klen_inc;
        uint32_t k2len = 2 * len - i + klen_inc;
        uint32_t k4len = 4 * len - i + klen_inc;
        for (uint32_t j = i; j < i + (i_inc / 8); j += 2) {
          // Load
          // First Quarter
          float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
          float64x2_t vj0_img = vld1q_f64(&v->img[j]);
          float64x2_t vk0_rel = vld1q_f64(&v->rel[j + len]);
          float64x2_t vk0_img = vld1q_f64(&v->img[j + len]);
          // Second Quarter
          float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 2 * len]);
          float64x2_t vj1_img = vld1q_f64(&v->img[j + 2 * len]);
          float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 3 * len]);
          float64x2_t vk1_img = vld1q_f64(&v->img[j + 3 * len]);
          // Third Quarter
          float64x2_t vj2_rel = vld1q_f64(&v->rel[j + 4 * len]);
          float64x2_t vj2_img = vld1q_f64(&v->img[j + 4 * len]);
          float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 5 * len]);
          float64x2_t vk2_img = vld1q_f64(&v->img[j + 5 * len]);
          // Forth Quarter
          float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 6 * len]);
          float64x2_t vj3_img = vld1q_f64(&v->img[j + 6 * len]);
          float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 7 * len]);
          float64x2_t vk3_img = vld1q_f64(&v->img[j + 7 * len]);
          float64x2_t vt_rel, vt_img;
          // First Layer: two roots
          float64x2_t vw1_rel = vld1q_f64(&w->rel[klen + j]);
          float64x2_t vw1_img = vld1q_f64(&w->img[klen + j]);
          CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw1_rel, vw1_img,
                       vt_rel, vt_img);
          CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw1_rel, vw1_img,
                       vt_rel, vt_img);
          CT_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw1_rel, vw1_img,
                       vt_rel, vt_img);
          CT_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw1_rel, vw1_img,
                       vt_rel, vt_img);
          // Second Layer: four roots
          float64x2_t vw20_rel = vld1q_f64(&w->rel[k2len + j]);
          float64x2_t vw20_img = vld1q_f64(&w->img[k2len + j]);
          CT_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw20_rel, vw20_img,
                       vt_rel, vt_img);
          CT_BUTTERFLY(vj2_rel, vj2_img, vj3_rel, vj3_img, vw20_rel, vw20_img,
                       vt_rel, vt_img);
          float64x2_t vw21_rel = vld1q_f64(&w->rel[k2len + j + len]);
          float64x2_t vw21_img = vld1q_f64(&w->img[k2len + j + len]);
          CT_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw21_rel, vw21_img,
                       vt_rel, vt_img);
          CT_BUTTERFLY(vk2_rel, vk2_img, vk3_rel, vk3_img, vw21_rel, vw21_img,
                       vt_rel, vt_img);
          // Third Layer: eight roots
          float64x2_t vw30_rel = vld1q_f64(&w->rel[k4len + j]);
          float64x2_t vw30_img = vld1q_f64(&w->img[k4len + j]);
          CT_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw30_rel, vw30_img,
                       vt_rel, vt_img);
          float64x2_t vw31_rel = vld1q_f64(&w->rel[k4len + j + len]);
          float64x2_t vw31_img = vld1q_f64(&w->img[k4len + j + len]);
          CT_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw31_rel, vw31_img,
                       vt_rel, vt_img);
          float64x2_t vw32_rel = vld1q_f64(&w->rel[k4len + j + 2 * len]);
          float64x2_t vw32_img = vld1q_f64(&w->img[k4len + j + 2 * len]);
          CT_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw32_rel, vw32_img,
                       vt_rel, vt_img);
          float64x2_t vw33_rel = vld1q_f64(&w->rel[k4len + j + 3 * len]);
          float64x2_t vw33_img = vld1q_f64(&w->img[k4len + j + 3 * len]);
          CT_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw33_rel, vw33_img,
                       vt_rel, vt_img);
          // Store
          // First Quarter
          vst1q_f64(&v->rel[j], vj0_rel);
          vst1q_f64(&v->img[j], vj0_img);
          vst1q_f64(&v->rel[j + len], vk0_rel);
          vst1q_f64(&v->img[j + len], vk0_img);
          // Second Quarter
          vst1q_f64(&v->rel[j + 2 * len], vj1_rel);
          vst1q_f64(&v->img[j + 2 * len], vj1_img);
          vst1q_f64(&v->rel[j + 3 * len], vk1_rel);
          vst1q_f64(&v->img[j + 3 * len], vk1_img);
          // Third Quarter
          vst1q_f64(&v->rel[j + 4 * len], vj2_rel);
          vst1q_f64(&v->img[j + 4 * len], vj2_img);
          vst1q_f64(&v->rel[j + 5 * len], vk2_rel);
          vst1q_f64(&v->img[j + 5 * len], vk2_img);
          // Forth Quarter
          vst1q_f64(&v->rel[j + 6 * len], vj3_rel);
          vst1q_f64(&v->img[j + 6 * len], vj3_img);
          vst1q_f64(&v->rel[j + 7 * len], vk3_rel);
          vst1q_f64(&v->img[j + 7 * len], vk3_img);
        }
      }
    }
  }

  if (m3_rem == 0 || sync_before_rem) {
#pragma omp barrier
  }

  if (m3_rem == 2) {
    // Merge two layers
#pragma omp for
    for (uint32_t j = 0; j < m / 4; j += 4) {
      // Load
      // First Quarter
      float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
      float64x2_t vj0_img = vld1q_f64(&v->img[j]);
      float64x2_t vk0_rel = vld1q_f64(&v->rel[j + m / 4]);
      float64x2_t vk0_img = vld1q_f64(&v->img[j + m / 4]);
      // Second Quarter
      float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 2]);
      float64x2_t vj1_img = vld1q_f64(&v->img[j + 2]);
      float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 2 + m / 4]);
      float64x2_t vk1_img = vld1q_f64(&v->img[j + 2 + m / 4]);
      // Third Quarter
      float64x2_t vj2_rel = vld1q_f64(&v->rel[j + m / 2]);
      float64x2_t vj2_img = vld1q_f64(&v->img[j + m / 2]);
      float64x2_t vk2_rel = vld1q_f64(&v->rel[j + m / 2 + m / 4]);
      float64x2_t vk2_img = vld1q_f64(&v->img[j + m / 2 + m / 4]);
      // Forth Quarter
      float64x2_t vj3_rel = vld1q_f64(&v->rel[j + m / 2 + 2]);
      float64x2_t vj3_img = vld1q_f64(&v->img[j + m / 2 + 2]);
      float64x2_t vk3_rel = vld1q_f64(&v->rel[j + m / 2 + 2 + m / 4]);
      float64x2_t vk3_img = vld1q_f64(&v->img[j + m / 2 + 2 + m / 4]);
      // First Layer
      float64x2_t vt_rel, vt_img;
      float64x2_t vw_rel, vw_img;
      vw_rel = vld1q_f64(&w->rel[m / 4 + j]);
      vw_img = vld1q_f64(&w->img[m / 4 + j]);
      CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      CT_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 4 + j + 2]);
      vw_img = vld1q_f64(&w->img[m / 4 + j + 2]);
      CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      CT_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Second Layer
      vw_rel = vld1q_f64(&w->rel[m / 2 + j]);
      vw_img = vld1q_f64(&w->img[m / 2 + j]);
      CT_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + 2]);
      vw_img = vld1q_f64(&w->img[m / 2 + j + 2]);
      CT_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + m / 4]);
      vw_img = vld1q_f64(&w->img[m / 2 + j + m / 4]);
      CT_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + m / 4 + 2]);
      vw_img = vld1q_f64(&w->img[m / 2 + j + m / 4 + 2]);
      CT_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Store
      // First Quarter
      vst1q_f64(&v->rel[j], vj0_rel);
      vst1q_f64(&v->img[j], vj0_img);
      vst1q_f64(&v->rel[j + m / 4], vk0_rel);
      vst1q_f64(&v->img[j + m / 4], vk0_img);
      // Second Quarter
      vst1q_f64(&v->rel[j + 2], vj1_rel);
      vst1q_f64(&v->img[j + 2], vj1_img);
      vst1q_f64(&v->rel[j + m / 4 + 2], vk1_rel);
      vst1q_f64(&v->img[j + m / 4 + 2], vk1_img);
      // Third Quarter
      vst1q_f64(&v->rel[j + m / 2], vj2_rel);
      vst1q_f64(&v->img[j + m / 2], vj2_img);
      vst1q_f64(&v->rel[j + m / 2 + m / 4], vk2_rel);
      vst1q_f64(&v->img[j + m / 2 + m / 4], vk2_img);
      // Forth Quarter
      vst1q_f64(&v->rel[j + m / 2 + 2], vj3_rel);
      vst1q_f64(&v->img[j + m / 2 + 2], vj3_img);
      vst1q_f64(&v->rel[j + m / 2 + m / 4 + 2], vk3_rel);
      vst1q_f64(&v->img[j + m / 2 + m / 4 + 2], vk3_img);
    }
  } else if (m3_rem == 1) {
    // Last layer
    if (chunk_size >= m / 2) chunk_size = (m / 2) / n_threads;
    for (uint32_t jj = tid * chunk_size; jj < m / 2;
         jj += n_threads * chunk_size) {
      for (uint32_t j = jj; j < jj + chunk_size; j += 8) {
        // Load
        // First Quarter
        float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
        float64x2_t vj0_img = vld1q_f64(&v->img[j]);
        float64x2_t vk0_rel = vld1q_f64(&v->rel[j + m / 2]);
        float64x2_t vk0_img = vld1q_f64(&v->img[j + m / 2]);
        // Second Quarter
        float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 2]);
        float64x2_t vj1_img = vld1q_f64(&v->img[j + 2]);
        float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 2 + m / 2]);
        float64x2_t vk1_img = vld1q_f64(&v->img[j + 2 + m / 2]);
        // Third Quarter
        float64x2_t vj2_rel = vld1q_f64(&v->rel[j + 4]);
        float64x2_t vj2_img = vld1q_f64(&v->img[j + 4]);
        float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 4 + m / 2]);
        float64x2_t vk2_img = vld1q_f64(&v->img[j + 4 + m / 2]);
        // Forth Quarter
        float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 6]);
        float64x2_t vj3_img = vld1q_f64(&v->img[j + 6]);
        float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 6 + m / 2]);
        float64x2_t vk3_img = vld1q_f64(&v->img[j + 6 + m / 2]);
        // CT
        float64x2_t vt_rel, vt_img;
        float64x2_t vw_rel, vw_img;
        vw_rel = vld1q_f64(&w->rel[m / 2 + j]);
        vw_img = vld1q_f64(&w->img[m / 2 + j]);
        CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 2]);
        vw_img = vld1q_f64(&w->img[m / 2 + j + 2]);
        CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 4]);
        vw_img = vld1q_f64(&w->img[m / 2 + j + 4]);
        CT_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 6]);
        vw_img = vld1q_f64(&w->img[m / 2 + j + 6]);
        CT_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        // Store
        // First Quarter
        vst1q_f64(&v->rel[j], vj0_rel);
        vst1q_f64(&v->img[j], vj0_img);
        vst1q_f64(&v->rel[j + m / 2], vk0_rel);
        vst1q_f64(&v->img[j + m / 2], vk0_img);
        // Second Quarter
        vst1q_f64(&v->rel[j + 2], vj1_rel);
        vst1q_f64(&v->img[j + 2], vj1_img);
        vst1q_f64(&v->rel[j + m / 2 + 2], vk1_rel);
        vst1q_f64(&v->img[j + m / 2 + 2], vk1_img);
        // Third Quarter
        vst1q_f64(&v->rel[j + 4], vj2_rel);
        vst1q_f64(&v->img[j + 4], vj2_img);
        vst1q_f64(&v->rel[j + m / 2 + 4], vk2_rel);
        vst1q_f64(&v->img[j + m / 2 + 4], vk2_img);
        // Forth Quarter
        vst1q_f64(&v->rel[j + 6], vj3_rel);
        vst1q_f64(&v->img[j + 6], vj3_img);
        vst1q_f64(&v->rel[j + m / 2 + 6], vk3_rel);
        vst1q_f64(&v->img[j + m / 2 + 6], vk3_img);
      }
    }
#pragma omp barrier
  }
#elif __AVX2_MODE__
  if (m < 16) {
#pragma omp single
    for (uint32_t len = 1; len < m; len <<= 1) {
      for (uint32_t i = 0; i < m; i += 2 * len) {
        for (uint32_t j = i; j < i + len; j++) {
          double w_rel = w->rel[len + j - i];
          double w_img = w->img[len + j - i];
          double j_rel = v->rel[j];
          double j_img = v->img[j];
          double k_rel = v->rel[j + len];
          double k_img = v->img[j + len];
          double t_rel = rel_mul(k_rel, k_img, w_rel, w_img);
          double t_img = img_mul(k_rel, k_img, w_rel, w_img);
          v->rel[j + len] = j_rel - t_rel;
          v->img[j + len] = j_img - t_img;
          v->rel[j] = j_rel + t_rel;
          v->img[j] = j_img + t_img;
        }
      }
    }
    return;
  }
  // Merge 4 layers using all the 16 AVX2 registers
  uint32_t k1 = 1;
  uint32_t k2 = 2;
  uint32_t k3 = 4;
  uint32_t k4 = 8;
  double w1_rel = w->rel[k1];
  double w1_img = w->img[k1];
  __m256d vw2_rel =
      _mm256_setr_pd(w->rel[k2], w->rel[k2 + 1], w->rel[k2], w->rel[k2 + 1]);
  __m256d vw2_img =
      _mm256_setr_pd(w->img[k2], w->img[k2 + 1], w->img[k2], w->img[k2 + 1]);
#pragma omp for
  for (uint32_t j = 0; j < m; j += 16) {
    // Load
    // First Half
    __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
    __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
    __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + 4]);
    __m256d vk0_img = _mm256_loadu_pd(&v->img[j + 4]);
    // Second Half
    __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 8]);
    __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 8]);
    __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 12]);
    __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 12]);
    // Tranpose
    // First Half
    __m256d vt0_rel = _mm256_permute2f128_pd(vj0_rel, vk0_rel, 0x20);
    __m256d vt0_img = _mm256_permute2f128_pd(vj0_img, vk0_img, 0x20);
    vk0_rel = _mm256_permute2f128_pd(vj0_rel, vk0_rel, 0x31);
    vk0_img = _mm256_permute2f128_pd(vj0_img, vk0_img, 0x31);
    vj0_rel = _mm256_unpacklo_pd(vt0_rel, vk0_rel);
    vj0_img = _mm256_unpacklo_pd(vt0_img, vk0_img);
    vk0_rel = _mm256_unpackhi_pd(vt0_rel, vk0_rel);
    vk0_img = _mm256_unpackhi_pd(vt0_img, vk0_img);
    // Second Half
    __m256d vt1_rel = _mm256_permute2f128_pd(vj1_rel, vk1_rel, 0x20);
    __m256d vt1_img = _mm256_permute2f128_pd(vj1_img, vk1_img, 0x20);
    vk1_rel = _mm256_permute2f128_pd(vj1_rel, vk1_rel, 0x31);
    vk1_img = _mm256_permute2f128_pd(vj1_img, vk1_img, 0x31);
    vj1_rel = _mm256_unpacklo_pd(vt1_rel, vk1_rel);
    vj1_img = _mm256_unpacklo_pd(vt1_img, vk1_img);
    vk1_rel = _mm256_unpackhi_pd(vt1_rel, vk1_rel);
    vk1_img = _mm256_unpackhi_pd(vt1_img, vk1_img);
    // First layer (len = 1), four roots per half
    // First Half
    __m256d vw_rel = _mm256_set1_pd(w1_rel);
    __m256d vw_img = _mm256_set1_pd(w1_img);
    CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Transpose back
    // First Half
    vt0_rel = _mm256_unpacklo_pd(vj0_rel, vk0_rel);
    vt0_img = _mm256_unpacklo_pd(vj0_img, vk0_img);
    vk0_rel = _mm256_unpackhi_pd(vj0_rel, vk0_rel);
    vk0_img = _mm256_unpackhi_pd(vj0_img, vk0_img);
    // Second Half
    vt1_rel = _mm256_unpacklo_pd(vj1_rel, vk1_rel);
    vt1_img = _mm256_unpacklo_pd(vj1_img, vk1_img);
    vk1_rel = _mm256_unpackhi_pd(vj1_rel, vk1_rel);
    vk1_img = _mm256_unpackhi_pd(vj1_img, vk1_img);
    // Second layer (len = 2), two roots per half
    // First Half
    CT_BUTTERFLY(vt0_rel, vt0_img, vk0_rel, vk0_img, vw2_rel, vw2_img, vj0_rel,
                 vj0_img);
    // Second Half
    CT_BUTTERFLY(vt1_rel, vt1_img, vk1_rel, vk1_img, vw2_rel, vw2_img, vj1_rel,
                 vj1_img);
    // Transpose back
    // First Half
    vj0_rel = _mm256_permute2f128_pd(vt0_rel, vk0_rel, 0x20);
    vj0_img = _mm256_permute2f128_pd(vt0_img, vk0_img, 0x20);
    vk0_rel = _mm256_permute2f128_pd(vt0_rel, vk0_rel, 0x31);
    vk0_img = _mm256_permute2f128_pd(vt0_img, vk0_img, 0x31);
    // Second Half
    vj1_rel = _mm256_permute2f128_pd(vt1_rel, vk1_rel, 0x20);
    vj1_img = _mm256_permute2f128_pd(vt1_img, vk1_img, 0x20);
    vk1_rel = _mm256_permute2f128_pd(vt1_rel, vk1_rel, 0x31);
    vk1_img = _mm256_permute2f128_pd(vt1_img, vk1_img, 0x31);
    // Third layer (len = 4), one root per half
    // First Half
    vw_rel = _mm256_loadu_pd(&w->rel[k3]);
    vw_img = _mm256_loadu_pd(&w->img[k3]);
    CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Forth layer (len = 8), one root
    // First Half
    vw_rel = _mm256_loadu_pd(&w->rel[k4]);
    vw_img = _mm256_loadu_pd(&w->img[k4]);
    CT_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    vw_rel = _mm256_loadu_pd(&w->rel[k4 + 4]);
    vw_img = _mm256_loadu_pd(&w->img[k4 + 4]);
    CT_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Store
    // First Half
    _mm256_storeu_pd(&v->rel[j], vj0_rel);
    _mm256_storeu_pd(&v->img[j], vj0_img);
    _mm256_storeu_pd(&v->rel[j + 4], vk0_rel);
    _mm256_storeu_pd(&v->img[j + 4], vk0_img);
    // Second Half
    _mm256_storeu_pd(&v->rel[j + 8], vj1_rel);
    _mm256_storeu_pd(&v->img[j + 8], vj1_img);
    _mm256_storeu_pd(&v->rel[j + 12], vk1_rel);
    _mm256_storeu_pd(&v->img[j + 12], vk1_img);
  }

  // Merge 2 consecutive layers
  uint32_t log2_m = floor_log2(m);
  uint32_t m2 = ((log2_m & 0x1) == 0) ? m : (m / 2);
  for (uint32_t len = 16; len < m2; len <<= 2) {
#pragma omp for
    for (uint32_t i = 0; i < m; i += 4 * len) {
      uint32_t klen = len - i;
      uint32_t k2len = 2 * len - i;
      for (uint32_t j = i; j < i + len; j += 4) {
        // Load
        // First Half
        __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
        __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
        __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + len]);
        __m256d vk0_img = _mm256_loadu_pd(&v->img[j + len]);
        // Second Half
        __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 2 * len]);
        __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 2 * len]);
        __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 3 * len]);
        __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 3 * len]);
        // We only have the upper half roots, so we conjugate here
        // First Layer: four roots
        __m256d vw1_rel = _mm256_loadu_pd(&w->rel[klen + j]);
        __m256d vw1_img = _mm256_loadu_pd(&w->img[klen + j]);
        // Second Layer: eight roots
        __m256d vw20_rel = _mm256_loadu_pd(&w->rel[k2len + j]);
        __m256d vw20_img = _mm256_loadu_pd(&w->img[k2len + j]);
        __m256d vw21_rel = _mm256_loadu_pd(&w->rel[k2len + j + len]);
        __m256d vw21_img = _mm256_loadu_pd(&w->img[k2len + j + len]);
        // First Layer (len)
        __m256d vt_rel, vt_img;
        CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw1_rel, vw1_img,
                     vt_rel, vt_img);
        CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw1_rel, vw1_img,
                     vt_rel, vt_img);
        // Second layer (2*len)
        CT_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw20_rel, vw20_img,
                     vt_rel, vt_img);
        CT_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw21_rel, vw21_img,
                     vt_rel, vt_img);
        // Store
        // First Half
        _mm256_storeu_pd(&v->rel[j], vj0_rel);
        _mm256_storeu_pd(&v->img[j], vj0_img);
        _mm256_storeu_pd(&v->rel[j + len], vk0_rel);
        _mm256_storeu_pd(&v->img[j + len], vk0_img);
        // Second Half
        _mm256_storeu_pd(&v->rel[j + 2 * len], vj1_rel);
        _mm256_storeu_pd(&v->img[j + 2 * len], vj1_img);
        _mm256_storeu_pd(&v->rel[j + 3 * len], vk1_rel);
        _mm256_storeu_pd(&v->img[j + 3 * len], vk1_img);
      }
    }
  }

  // Final layer for odd log2_n
  if (m2 < m) {
    // Unroll two vector iterations
#pragma omp for
    for (uint32_t j = 0; j < m / 2; j += 8) {
      // Load
      // First Half
      __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
      __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
      __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + m / 2]);
      __m256d vk0_img = _mm256_loadu_pd(&v->img[j + m / 2]);
      // Second Half
      __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 4]);
      __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 4]);
      __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 4 + m / 2]);
      __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 4 + m / 2]);
      // CT
      __m256d vt_rel, vt_img;
      __m256d vw_rel, vw_img;
      vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j]);
      vw_img = _mm256_loadu_pd(&w->img[m / 2 + j]);
      CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j + 4]);
      vw_img = _mm256_loadu_pd(&w->img[m / 2 + j + 4]);
      CT_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Store
      // First Half
      _mm256_storeu_pd(&v->rel[j], vj0_rel);
      _mm256_storeu_pd(&v->img[j], vj0_img);
      _mm256_storeu_pd(&v->rel[j + m / 2], vk0_rel);
      _mm256_storeu_pd(&v->img[j + m / 2], vk0_img);
      // Second Half
      _mm256_storeu_pd(&v->rel[j + 4], vj1_rel);
      _mm256_storeu_pd(&v->img[j + 4], vj1_img);
      _mm256_storeu_pd(&v->rel[j + 4 + m / 2], vk1_rel);
      _mm256_storeu_pd(&v->img[j + 4 + m / 2], vk1_img);
    }
  }
#else
  for (uint32_t len = 1; len < m; len <<= 1) {
#pragma omp for
    for (uint32_t i = 0; i < m; i += 2 * len) {
      for (uint32_t j = i; j < i + len; j++) {
        double w_rel = w->rel[len + j - i];
        double w_img = w->img[len + j - i];
        double j_rel = v->rel[j];
        double j_img = v->img[j];
        double k_rel = v->rel[j + len];
        double k_img = v->img[j + len];
        double t_rel = rel_mul(k_rel, k_img, w_rel, w_img);
        double t_img = img_mul(k_rel, k_img, w_rel, w_img);
        v->rel[j + len] = j_rel - t_rel;
        v->img[j + len] = j_img - t_img;
        v->rel[j] = j_rel + t_rel;
        v->img[j] = j_img + t_img;
      }
    }
  }
#endif
}

#define GS_BUTTERFLY(vj_rel, vj_img, vk_rel, vk_img, vw_rel, vw_img, vt_rel, \
                     vt_img)                                                 \
  do {                                                                       \
    (vt_rel) = complex_vsub((vj_rel), (vk_rel));                             \
    (vt_img) = complex_vsub((vj_img), (vk_img));                             \
    (vj_rel) = complex_vadd((vj_rel), (vk_rel));                             \
    (vj_img) = complex_vadd((vj_img), (vk_img));                             \
    (vk_rel) = rel_vmul((vt_rel), (vt_img), (vw_rel), (vw_img));             \
    (vk_img) = img_vmul((vt_rel), (vt_img), (vw_rel), (vw_img));             \
  } while (0)

#ifdef __NEON_MODE__
#define GS_BUTTERFLY_X(vj_rel, vj_img, vk_rel, vk_img, w_rel, w_img, vt_rel, \
                       vt_img)                                               \
  do {                                                                       \
    (vt_rel) = complex_vsub((vj_rel), (vk_rel));                             \
    (vt_img) = complex_vsub((vj_img), (vk_img));                             \
    (vj_rel) = complex_vadd((vj_rel), (vk_rel));                             \
    (vj_img) = complex_vadd((vj_img), (vk_img));                             \
    (vk_rel) = rel_vmulx((vt_rel), (vt_img), (w_rel), (w_img));              \
    (vk_img) = img_vmulx((vt_rel), (vt_img), (w_rel), (w_img));              \
  } while (0)
#endif

// m = n/2 (half cyclotomic size)
void ifft(vec_t* v, const vec_t* w, uint32_t m) {
#ifdef __NEON_MODE__

  if (m < 16) {
#pragma omp single
    for (uint32_t len = m / 2; len > 0; len >>= 1) {
      for (uint32_t i = 0; i < m; i += 2 * len) {
        for (uint32_t j = i; j < i + len; j++) {
          double w_rel = w->rel[len + j - i];
          double w_img = -w->img[len + j - i];  // Conjugate
          double j_rel = v->rel[j];
          double j_img = v->img[j];
          double k_rel = v->rel[j + len];
          double k_img = v->img[j + len];
          double t_rel = j_rel - k_rel;
          double t_img = j_img - k_img;
          v->rel[j] = j_rel + k_rel;
          v->img[j] = j_img + k_img;
          v->rel[j + len] = rel_mul(t_rel, t_img, w_rel, w_img);
          v->img[j + len] = img_mul(t_rel, t_img, w_rel, w_img);
        }
      }
    }
    return;
  }

  // Based on the number of threads and the format of n,
  // we decide the almost-optimal block-division
  uint32_t n_threads = (uint32_t)omp_get_num_threads();
  uint32_t tid = (uint32_t)omp_get_thread_num();
  uint32_t chunk_size = m / n_threads;
  uint32_t log2_m = floor_log2(m / 16);
  uint32_t m3_rem = log2_m % 3;
  uint32_t m3 = 1 << (log2_m + 4 - m3_rem);
  bool sync_first_3x = chunk_size < m3 && chunk_size < m;
  if (sync_first_3x) chunk_size = m3 >> 3;

  if (m3_rem == 2) {
    // Merge two layers
#pragma omp for
    for (uint32_t j = 0; j < m / 4; j += 4) {
      // Load
      // First Quarter
      float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
      float64x2_t vj0_img = vld1q_f64(&v->img[j]);
      float64x2_t vk0_rel = vld1q_f64(&v->rel[j + m / 4]);
      float64x2_t vk0_img = vld1q_f64(&v->img[j + m / 4]);
      // Second Quarter
      float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 2]);
      float64x2_t vj1_img = vld1q_f64(&v->img[j + 2]);
      float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 2 + m / 4]);
      float64x2_t vk1_img = vld1q_f64(&v->img[j + 2 + m / 4]);
      // Third Quarter
      float64x2_t vj2_rel = vld1q_f64(&v->rel[j + m / 2]);
      float64x2_t vj2_img = vld1q_f64(&v->img[j + m / 2]);
      float64x2_t vk2_rel = vld1q_f64(&v->rel[j + m / 2 + m / 4]);
      float64x2_t vk2_img = vld1q_f64(&v->img[j + m / 2 + m / 4]);
      // Forth Quarter
      float64x2_t vj3_rel = vld1q_f64(&v->rel[j + m / 2 + 2]);
      float64x2_t vj3_img = vld1q_f64(&v->img[j + m / 2 + 2]);
      float64x2_t vk3_rel = vld1q_f64(&v->rel[j + m / 2 + 2 + m / 4]);
      float64x2_t vk3_img = vld1q_f64(&v->img[j + m / 2 + 2 + m / 4]);
      // First Layer
      float64x2_t vt_rel, vt_img;
      float64x2_t vw_rel, vw_img;
      vw_rel = vld1q_f64(&w->rel[m / 2 + j]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j]));
      GS_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + 2]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + 2]));
      GS_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + m / 4]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + m / 4]));
      GS_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 2 + j + m / 4 + 2]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + m / 4 + 2]));
      GS_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Second Layer
      vw_rel = vld1q_f64(&w->rel[m / 4 + j]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 4 + j]));
      GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      GS_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = vld1q_f64(&w->rel[m / 4 + j + 2]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[m / 4 + j + 2]));
      GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      GS_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Store
      // First Quarter
      vst1q_f64(&v->rel[j], vj0_rel);
      vst1q_f64(&v->img[j], vj0_img);
      vst1q_f64(&v->rel[j + m / 4], vk0_rel);
      vst1q_f64(&v->img[j + m / 4], vk0_img);
      // Second Quarter
      vst1q_f64(&v->rel[j + 2], vj1_rel);
      vst1q_f64(&v->img[j + 2], vj1_img);
      vst1q_f64(&v->rel[j + m / 4 + 2], vk1_rel);
      vst1q_f64(&v->img[j + m / 4 + 2], vk1_img);
      // Third Quarter
      vst1q_f64(&v->rel[j + m / 2], vj2_rel);
      vst1q_f64(&v->img[j + m / 2], vj2_img);
      vst1q_f64(&v->rel[j + m / 2 + m / 4], vk2_rel);
      vst1q_f64(&v->img[j + m / 2 + m / 4], vk2_img);
      // Forth Quarter
      vst1q_f64(&v->rel[j + m / 2 + 2], vj3_rel);
      vst1q_f64(&v->img[j + m / 2 + 2], vj3_img);
      vst1q_f64(&v->rel[j + m / 2 + m / 4 + 2], vk3_rel);
      vst1q_f64(&v->img[j + m / 2 + m / 4 + 2], vk3_img);
    }
  } else if (m3_rem == 1) {
    // Last layer
    uint32_t chunk_size2 = chunk_size;
    if (chunk_size2 >= m / 2) chunk_size2 = (m / 2) / n_threads;
    for (uint32_t jj = tid * chunk_size2; jj < m / 2;
         jj += n_threads * chunk_size2) {
      for (uint32_t j = jj; j < jj + chunk_size2; j += 8) {
        // Load
        // First Quarter
        float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
        float64x2_t vj0_img = vld1q_f64(&v->img[j]);
        float64x2_t vk0_rel = vld1q_f64(&v->rel[j + m / 2]);
        float64x2_t vk0_img = vld1q_f64(&v->img[j + m / 2]);
        // Second Quarter
        float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 2]);
        float64x2_t vj1_img = vld1q_f64(&v->img[j + 2]);
        float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 2 + m / 2]);
        float64x2_t vk1_img = vld1q_f64(&v->img[j + 2 + m / 2]);
        // Third Quarter
        float64x2_t vj2_rel = vld1q_f64(&v->rel[j + 4]);
        float64x2_t vj2_img = vld1q_f64(&v->img[j + 4]);
        float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 4 + m / 2]);
        float64x2_t vk2_img = vld1q_f64(&v->img[j + 4 + m / 2]);
        // Forth Quarter
        float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 6]);
        float64x2_t vj3_img = vld1q_f64(&v->img[j + 6]);
        float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 6 + m / 2]);
        float64x2_t vk3_img = vld1q_f64(&v->img[j + 6 + m / 2]);
        // CT
        float64x2_t vt_rel, vt_img;
        float64x2_t vw_rel, vw_img;
        vw_rel = vld1q_f64(&w->rel[m / 2 + j]);
        vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j]));
        GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 2]);
        vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + 2]));
        GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 4]);
        vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + 4]));
        GS_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        vw_rel = vld1q_f64(&w->rel[m / 2 + j + 6]);
        vw_img = vnegq_f64(vld1q_f64(&w->img[m / 2 + j + 6]));
        GS_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw_rel, vw_img, vt_rel,
                     vt_img);
        // Store
        // First Quarter
        vst1q_f64(&v->rel[j], vj0_rel);
        vst1q_f64(&v->img[j], vj0_img);
        vst1q_f64(&v->rel[j + m / 2], vk0_rel);
        vst1q_f64(&v->img[j + m / 2], vk0_img);
        // Second Quarter
        vst1q_f64(&v->rel[j + 2], vj1_rel);
        vst1q_f64(&v->img[j + 2], vj1_img);
        vst1q_f64(&v->rel[j + m / 2 + 2], vk1_rel);
        vst1q_f64(&v->img[j + m / 2 + 2], vk1_img);
        // Third Quarter
        vst1q_f64(&v->rel[j + 4], vj2_rel);
        vst1q_f64(&v->img[j + 4], vj2_img);
        vst1q_f64(&v->rel[j + m / 2 + 4], vk2_rel);
        vst1q_f64(&v->img[j + m / 2 + 4], vk2_img);
        // Forth Quarter
        vst1q_f64(&v->rel[j + 6], vj3_rel);
        vst1q_f64(&v->img[j + 6], vj3_img);
        vst1q_f64(&v->rel[j + m / 2 + 6], vk3_rel);
        vst1q_f64(&v->img[j + m / 2 + 6], vk3_img);
      }
    }
#pragma omp barrier
  }

  // Merge 3 consecutive layers
  uint32_t jj_inc = n_threads * chunk_size;
  for (uint32_t len = m3 / 2; len > 8; len >>= 3) {
    uint32_t i_inc = 2 * len;
    if (sync_first_3x && len == (m3 >> 1)) {
      jj_inc = m3;
      i_inc = 8 * chunk_size;
    }
    for (uint32_t jj = tid * chunk_size; jj < m; jj += jj_inc) {
      for (uint32_t i = 0; i < m; i += 2 * len) {
        uint32_t klen = len - i;
        uint32_t k2len = len / 2 - i;
        uint32_t k4len = len / 4 - i;
        for (uint32_t j = i; j < i + (i_inc / 8); j += 2) {
          // Load
          // First Quarter
          float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
          float64x2_t vj0_img = vld1q_f64(&v->img[j]);
          float64x2_t vk0_rel = vld1q_f64(&v->rel[j + len]);
          float64x2_t vk0_img = vld1q_f64(&v->img[j + len]);
          // Second Quarter
          float64x2_t vj1_rel = vld1q_f64(&v->rel[j + len / 2]);
          float64x2_t vj1_img = vld1q_f64(&v->img[j + len / 2]);
          float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 3 * len / 2]);
          float64x2_t vk1_img = vld1q_f64(&v->img[j + 3 * len / 2]);
          // Third Quarter
          float64x2_t vj2_rel = vld1q_f64(&v->rel[j + len / 4]);
          float64x2_t vj2_img = vld1q_f64(&v->img[j + len / 4]);
          float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 5 * len / 4]);
          float64x2_t vk2_img = vld1q_f64(&v->img[j + 5 * len / 4]);
          // Forth Quarter
          float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 3 * len / 4]);
          float64x2_t vj3_img = vld1q_f64(&v->img[j + 3 * len / 4]);
          float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 7 * len / 4]);
          float64x2_t vk3_img = vld1q_f64(&v->img[j + 7 * len / 4]);
          float64x2_t vt_rel, vt_img;
          // First Layer: eight roots
          float64x2_t vw10_rel = vld1q_f64(&w->rel[klen + j]);
          float64x2_t vw10_img = vnegq_f64(vld1q_f64(&w->img[klen + j]));
          GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw10_rel, vw10_img,
                       vt_rel, vt_img);
          float64x2_t vw12_rel = vld1q_f64(&w->rel[klen + j + len / 4]);
          float64x2_t vw12_img =
              vnegq_f64(vld1q_f64(&w->img[klen + j + len / 4]));
          GS_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw12_rel, vw12_img,
                       vt_rel, vt_img);
          float64x2_t vw11_rel = vld1q_f64(&w->rel[klen + j + 2 * len / 4]);
          float64x2_t vw11_img =
              vnegq_f64(vld1q_f64(&w->img[klen + j + 2 * len / 4]));
          GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw11_rel, vw11_img,
                       vt_rel, vt_img);
          float64x2_t vw13_rel = vld1q_f64(&w->rel[klen + j + 3 * len / 4]);
          float64x2_t vw13_img =
              vnegq_f64(vld1q_f64(&w->img[klen + j + 3 * len / 4]));
          GS_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw13_rel, vw13_img,
                       vt_rel, vt_img);
          // Second Layer: four roots
          float64x2_t vw20_rel = vld1q_f64(&w->rel[k2len + j]);
          float64x2_t vw20_img = vnegq_f64(vld1q_f64(&w->img[k2len + j]));
          GS_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw20_rel, vw20_img,
                       vt_rel, vt_img);
          GS_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw20_rel, vw20_img,
                       vt_rel, vt_img);
          float64x2_t vw21_rel = vld1q_f64(&w->rel[k2len + j + len / 4]);
          float64x2_t vw21_img =
              vnegq_f64(vld1q_f64(&w->img[k2len + j + len / 4]));
          GS_BUTTERFLY(vj2_rel, vj2_img, vj3_rel, vj3_img, vw21_rel, vw21_img,
                       vt_rel, vt_img);
          GS_BUTTERFLY(vk2_rel, vk2_img, vk3_rel, vk3_img, vw21_rel, vw21_img,
                       vt_rel, vt_img);
          // Third Layer: two roots
          float64x2_t vw3_rel = vld1q_f64(&w->rel[k4len + j]);
          float64x2_t vw3_img = vnegq_f64(vld1q_f64(&w->img[k4len + j]));
          GS_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw3_rel, vw3_img,
                       vt_rel, vt_img);
          GS_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw3_rel, vw3_img,
                       vt_rel, vt_img);
          GS_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw3_rel, vw3_img,
                       vt_rel, vt_img);
          GS_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw3_rel, vw3_img,
                       vt_rel, vt_img);
          // Store
          // First Quarter
          vst1q_f64(&v->rel[j], vj0_rel);
          vst1q_f64(&v->img[j], vj0_img);
          vst1q_f64(&v->rel[j + len], vk0_rel);
          vst1q_f64(&v->img[j + len], vk0_img);
          // Second Quarter
          vst1q_f64(&v->rel[j + len / 2], vj1_rel);
          vst1q_f64(&v->img[j + len / 2], vj1_img);
          vst1q_f64(&v->rel[j + 3 * len / 2], vk1_rel);
          vst1q_f64(&v->img[j + 3 * len / 2], vk1_img);
          // Third Quarter
          vst1q_f64(&v->rel[j + len / 4], vj2_rel);
          vst1q_f64(&v->img[j + len / 4], vj2_img);
          vst1q_f64(&v->rel[j + 5 * len / 4], vk2_rel);
          vst1q_f64(&v->img[j + 5 * len / 4], vk2_img);
          // Forth Quarter
          vst1q_f64(&v->rel[j + 3 * len / 4], vj3_rel);
          vst1q_f64(&v->img[j + 3 * len / 4], vj3_img);
          vst1q_f64(&v->rel[j + 7 * len / 4], vk3_rel);
          vst1q_f64(&v->img[j + 7 * len / 4], vk3_img);
        }
        if (sync_first_3x && len == (m3 >> 1)) {
          jj_inc = n_threads * chunk_size;
#pragma omp barrier
        }
      }
    }
  }

  // Merge 4 layers using all the 32 ARM NEON (ARMv8) registers
  uint32_t k1 = 1;
  uint32_t k2 = 2;
  uint32_t k3 = 4;
  uint32_t k4 = 8;
  double w1_rel = w->rel[k1];
  double w1_img = -w->img[k1];
  float64x2_t vw2_rel = vld1q_f64(&w->rel[k2]);
  float64x2_t vw2_img = vnegq_f64(vld1q_f64(&w->img[k2]));
  for (uint32_t jj = tid * chunk_size; jj < m; jj += n_threads * chunk_size) {
    for (uint32_t j = jj; j < jj + chunk_size; j += 16) {
      // Load
      // First Quarter
      float64x2_t vj0_rel = vld1q_f64(&v->rel[j]);
      float64x2_t vj0_img = vld1q_f64(&v->img[j]);
      float64x2_t vk0_rel = vld1q_f64(&v->rel[j + 2]);
      float64x2_t vk0_img = vld1q_f64(&v->img[j + 2]);
      // Second Quarter
      float64x2_t vj1_rel = vld1q_f64(&v->rel[j + 4]);
      float64x2_t vj1_img = vld1q_f64(&v->img[j + 4]);
      float64x2_t vk1_rel = vld1q_f64(&v->rel[j + 6]);
      float64x2_t vk1_img = vld1q_f64(&v->img[j + 6]);
      // Third Quarter
      float64x2_t vj2_rel = vld1q_f64(&v->rel[j + 8]);
      float64x2_t vj2_img = vld1q_f64(&v->img[j + 8]);
      float64x2_t vk2_rel = vld1q_f64(&v->rel[j + 10]);
      float64x2_t vk2_img = vld1q_f64(&v->img[j + 10]);
      // Forth Quarter
      float64x2_t vj3_rel = vld1q_f64(&v->rel[j + 12]);
      float64x2_t vj3_img = vld1q_f64(&v->img[j + 12]);
      float64x2_t vk3_rel = vld1q_f64(&v->rel[j + 14]);
      float64x2_t vk3_img = vld1q_f64(&v->img[j + 14]);
      // First layer (len = 8), eight roots
      float64x2_t vw_rel, vw_img;
      // First Quarter
      float64x2_t vt0_rel, vt0_img;
      vw_rel = vld1q_f64(&w->rel[k4]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k4]));
      GS_BUTTERFLY(vj0_rel, vj0_img, vj2_rel, vj2_img, vw_rel, vw_img, vt0_rel,
                   vt0_img);
      // Second Quarter
      float64x2_t vt1_rel, vt1_img;
      vw_rel = vld1q_f64(&w->rel[k4 + 2]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k4 + 2]));
      GS_BUTTERFLY(vk0_rel, vk0_img, vk2_rel, vk2_img, vw_rel, vw_img, vt1_rel,
                   vt1_img);
      // Third Quarter
      float64x2_t vt2_rel, vt2_img;
      vw_rel = vld1q_f64(&w->rel[k4 + 4]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k4 + 4]));
      GS_BUTTERFLY(vj1_rel, vj1_img, vj3_rel, vj3_img, vw_rel, vw_img, vt2_rel,
                   vt2_img);
      // Quarter Quarter
      float64x2_t vt3_rel, vt3_img;
      vw_rel = vld1q_f64(&w->rel[k4 + 6]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k4 + 6]));
      GS_BUTTERFLY(vk1_rel, vk1_img, vk3_rel, vk3_img, vw_rel, vw_img, vt3_rel,
                   vt3_img);
      // Second layer (len = 4), four roots
      // First Quarter
      vw_rel = vld1q_f64(&w->rel[k3]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k3]));
      GS_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw_rel, vw_img, vt0_rel,
                   vt0_img);
      // Third Quarter
      GS_BUTTERFLY(vj2_rel, vj2_img, vj3_rel, vj3_img, vw_rel, vw_img, vt2_rel,
                   vt2_img);
      // Second Quarter
      vw_rel = vld1q_f64(&w->rel[k3 + 2]);
      vw_img = vnegq_f64(vld1q_f64(&w->img[k3 + 2]));
      GS_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                   vt1_img);
      // Forth Quarter
      GS_BUTTERFLY(vk2_rel, vk2_img, vk3_rel, vk3_img, vw_rel, vw_img, vt3_rel,
                   vt3_img);
      // Third layer (len = 2), two roots
      // First Quarter
      GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw2_rel, vw2_img,
                   vt0_rel, vt0_img);
      // Second Quarter
      GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw2_rel, vw2_img,
                   vt1_rel, vt1_img);
      // Third Quarter
      GS_BUTTERFLY(vj2_rel, vj2_img, vk2_rel, vk2_img, vw2_rel, vw2_img,
                   vt2_rel, vt2_img);
      // Forth Quarter
      GS_BUTTERFLY(vj3_rel, vj3_img, vk3_rel, vk3_img, vw2_rel, vw2_img,
                   vt3_rel, vt3_img);
      // Forth layer (len = 1), one root
      // Tranpose
      // First Quarter
      vt0_rel = vzip1q_f64(vj0_rel, vk0_rel);
      vt0_img = vzip1q_f64(vj0_img, vk0_img);
      vk0_rel = vzip2q_f64(vj0_rel, vk0_rel);
      vk0_img = vzip2q_f64(vj0_img, vk0_img);
      // Second Quarter
      vt1_rel = vzip1q_f64(vj1_rel, vk1_rel);
      vt1_img = vzip1q_f64(vj1_img, vk1_img);
      vk1_rel = vzip2q_f64(vj1_rel, vk1_rel);
      vk1_img = vzip2q_f64(vj1_img, vk1_img);
      // Third Quarter
      vt2_rel = vzip1q_f64(vj2_rel, vk2_rel);
      vt2_img = vzip1q_f64(vj2_img, vk2_img);
      vk2_rel = vzip2q_f64(vj2_rel, vk2_rel);
      vk2_img = vzip2q_f64(vj2_img, vk2_img);
      // Forth Quarter
      vt3_rel = vzip1q_f64(vj3_rel, vk3_rel);
      vt3_img = vzip1q_f64(vj3_img, vk3_img);
      vk3_rel = vzip2q_f64(vj3_rel, vk3_rel);
      vk3_img = vzip2q_f64(vj3_img, vk3_img);
      // First Quarter
      GS_BUTTERFLY_X(vt0_rel, vt0_img, vk0_rel, vk0_img, w1_rel, w1_img,
                     vj0_rel, vj0_img);
      // Second Quarter
      GS_BUTTERFLY_X(vt1_rel, vt1_img, vk1_rel, vk1_img, w1_rel, w1_img,
                     vj1_rel, vj1_img);
      // Third Quarter
      GS_BUTTERFLY_X(vt2_rel, vt2_img, vk2_rel, vk2_img, w1_rel, w1_img,
                     vj2_rel, vj2_img);
      // Forth Quarter
      GS_BUTTERFLY_X(vt3_rel, vt3_img, vk3_rel, vk3_img, w1_rel, w1_img,
                     vj3_rel, vj3_img);
      // Transpose back
      // First Quarter
      vj0_rel = vzip1q_f64(vt0_rel, vk0_rel);
      vj0_img = vzip1q_f64(vt0_img, vk0_img);
      vk0_rel = vzip2q_f64(vt0_rel, vk0_rel);
      vk0_img = vzip2q_f64(vt0_img, vk0_img);
      // Second Quarter
      vj1_rel = vzip1q_f64(vt1_rel, vk1_rel);
      vj1_img = vzip1q_f64(vt1_img, vk1_img);
      vk1_rel = vzip2q_f64(vt1_rel, vk1_rel);
      vk1_img = vzip2q_f64(vt1_img, vk1_img);
      // Third Quarter
      vj2_rel = vzip1q_f64(vt2_rel, vk2_rel);
      vj2_img = vzip1q_f64(vt2_img, vk2_img);
      vk2_rel = vzip2q_f64(vt2_rel, vk2_rel);
      vk2_img = vzip2q_f64(vt2_img, vk2_img);
      // Forth Quarter
      vj3_rel = vzip1q_f64(vt3_rel, vk3_rel);
      vj3_img = vzip1q_f64(vt3_img, vk3_img);
      vk3_rel = vzip2q_f64(vt3_rel, vk3_rel);
      vk3_img = vzip2q_f64(vt3_img, vk3_img);
      // Store
      // First Quarter
      vst1q_f64(&v->rel[j], vj0_rel);
      vst1q_f64(&v->img[j], vj0_img);
      vst1q_f64(&v->rel[j + 2], vk0_rel);
      vst1q_f64(&v->img[j + 2], vk0_img);
      // Second Quarter
      vst1q_f64(&v->rel[j + 4], vj1_rel);
      vst1q_f64(&v->img[j + 4], vj1_img);
      vst1q_f64(&v->rel[j + 6], vk1_rel);
      vst1q_f64(&v->img[j + 6], vk1_img);
      // Third Quarter
      vst1q_f64(&v->rel[j + 8], vj2_rel);
      vst1q_f64(&v->img[j + 8], vj2_img);
      vst1q_f64(&v->rel[j + 10], vk2_rel);
      vst1q_f64(&v->img[j + 10], vk2_img);
      // Forth Quarter
      vst1q_f64(&v->rel[j + 12], vj3_rel);
      vst1q_f64(&v->img[j + 12], vj3_img);
      vst1q_f64(&v->rel[j + 14], vk3_rel);
      vst1q_f64(&v->img[j + 14], vk3_img);
    }
  }

#pragma omp barrier

#elif __AVX2_MODE__
  if (m < 16) {
#pragma omp single
    for (uint32_t len = m / 2; len > 0; len >>= 1) {
      for (uint32_t i = 0; i < m; i += 2 * len) {
        for (uint32_t j = i; j < i + len; j++) {
          double w_rel = w->rel[len + j - i];
          double w_img = -w->img[len + j - i];  // Conjugate
          double j_rel = v->rel[j];
          double j_img = v->img[j];
          double k_rel = v->rel[j + len];
          double k_img = v->img[j + len];
          double t_rel = j_rel - k_rel;
          double t_img = j_img - k_img;
          v->rel[j] = j_rel + k_rel;
          v->img[j] = j_img + k_img;
          v->rel[j + len] = rel_mul(t_rel, t_img, w_rel, w_img);
          v->img[j + len] = img_mul(t_rel, t_img, w_rel, w_img);
        }
      }
    }
    return;
  }

  uint32_t log2_m = floor_log2(m);
  uint32_t m2 = ((log2_m & 0x1) == 0) ? m : (m / 2);
  __m256d sign = _mm256_set1_pd(-0.0f);
  // First layer for odd log2_m
  if (m2 < m) {
    // Unroll two vector iterations
#pragma omp for
    for (uint32_t j = 0; j < m / 2; j += 8) {
      // Load
      // First Half
      __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
      __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
      __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + m / 2]);
      __m256d vk0_img = _mm256_loadu_pd(&v->img[j + m / 2]);
      // Second Half
      __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 4]);
      __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 4]);
      __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 4 + m / 2]);
      __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 4 + m / 2]);
      // GS
      __m256d vt_rel, vt_img;
      __m256d vw_rel, vw_img;
      vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j]);
      vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[m / 2 + j]), sign);
      GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j + 4]);
      vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[m / 2 + j + 4]), sign);
      GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt_rel,
                   vt_img);
      // Store
      // First Half
      _mm256_storeu_pd(&v->rel[j], vj0_rel);
      _mm256_storeu_pd(&v->img[j], vj0_img);
      _mm256_storeu_pd(&v->rel[j + m / 2], vk0_rel);
      _mm256_storeu_pd(&v->img[j + m / 2], vk0_img);
      // Second Half
      _mm256_storeu_pd(&v->rel[j + 4], vj1_rel);
      _mm256_storeu_pd(&v->img[j + 4], vj1_img);
      _mm256_storeu_pd(&v->rel[j + 4 + m / 2], vk1_rel);
      _mm256_storeu_pd(&v->img[j + 4 + m / 2], vk1_img);
    }
  }

  // Merge 2 consecutive layers
  for (uint32_t len = m2 / 2; len > 8; len >>= 2) {
#pragma omp for
    for (uint32_t i = 0; i < m; i += 2 * len) {
      uint32_t klen = len - i;
      uint32_t k2len = len / 2 - i;
      for (uint32_t j = i; j < i + len / 2; j += 4) {
        // Load
        // First Half
        __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
        __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
        __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + len]);
        __m256d vk0_img = _mm256_loadu_pd(&v->img[j + len]);
        // Second Half
        __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + len / 2]);
        __m256d vj1_img = _mm256_loadu_pd(&v->img[j + len / 2]);
        __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 3 * len / 2]);
        __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 3 * len / 2]);
        // We only have the upper half roots, so we conjugate here
        // First Layer: eight roots
        __m256d vw10_rel = _mm256_loadu_pd(&w->rel[klen + j]);
        __m256d vw10_img =
            _mm256_xor_pd(_mm256_loadu_pd(&w->img[klen + j]), sign);
        __m256d vw11_rel = _mm256_loadu_pd(&w->rel[klen + j + len / 2]);
        __m256d vw11_img =
            _mm256_xor_pd(_mm256_loadu_pd(&w->img[klen + j + len / 2]), sign);
        // Second Layer: four roots
        __m256d vw2_rel = _mm256_loadu_pd(&w->rel[k2len + j]);
        __m256d vw2_img =
            _mm256_xor_pd(_mm256_loadu_pd(&w->img[k2len + j]), sign);
        // First Layer (len)
        __m256d vt_rel, vt_img;
        GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw10_rel, vw10_img,
                     vt_rel, vt_img);
        GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw11_rel, vw11_img,
                     vt_rel, vt_img);
        // Second layer (len/2)
        GS_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw2_rel, vw2_img,
                     vt_rel, vt_img);
        GS_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw2_rel, vw2_img,
                     vt_rel, vt_img);
        // Store
        // First Half
        _mm256_storeu_pd(&v->rel[j], vj0_rel);
        _mm256_storeu_pd(&v->img[j], vj0_img);
        _mm256_storeu_pd(&v->rel[j + len], vk0_rel);
        _mm256_storeu_pd(&v->img[j + len], vk0_img);
        // Second Half
        _mm256_storeu_pd(&v->rel[j + len / 2], vj1_rel);
        _mm256_storeu_pd(&v->img[j + len / 2], vj1_img);
        _mm256_storeu_pd(&v->rel[j + 3 * len / 2], vk1_rel);
        _mm256_storeu_pd(&v->img[j + 3 * len / 2], vk1_img);
      }
    }
  }

  // Merge 4 layers using all the 16 AVX2 registers
  uint32_t k1 = 1;
  uint32_t k2 = 2;
  uint32_t k3 = 4;
  uint32_t k4 = 8;
  double w1_rel = w->rel[k1];
  double w1_img = -w->img[k1];
  __m256d vw2_rel =
      _mm256_setr_pd(w->rel[k2], w->rel[k2 + 1], w->rel[k2], w->rel[k2 + 1]);
  __m256d vw2_img = _mm256_setr_pd(-w->img[k2], -w->img[k2 + 1], -w->img[k2],
                                   -w->img[k2 + 1]);
#pragma omp for
  for (uint32_t j = 0; j < m; j += 16) {
    // Load
    // First Half
    __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
    __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
    __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + 4]);
    __m256d vk0_img = _mm256_loadu_pd(&v->img[j + 4]);
    // Second Half
    __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 8]);
    __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 8]);
    __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 12]);
    __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 12]);
    // First layer (len = 8), one root
    __m256d vt0_rel, vt0_img;
    __m256d vt1_rel, vt1_img;
    __m256d vw_rel, vw_img;
    // First Half
    vw_rel = _mm256_loadu_pd(&w->rel[k4]);
    vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[k4]), sign);
    GS_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    vw_rel = _mm256_loadu_pd(&w->rel[k4 + 4]);
    vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[k4 + 4]), sign);
    GS_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Second layer (len = 4), one root per half
    // First Half
    vw_rel = _mm256_loadu_pd(&w->rel[k3]);
    vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[k3]), sign);
    GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Third layer (len = 2), two roots per half
    // Transpose
    // First Half
    vt0_rel = _mm256_permute2f128_pd(vj0_rel, vk0_rel, 0x20);
    vt0_img = _mm256_permute2f128_pd(vj0_img, vk0_img, 0x20);
    vk0_rel = _mm256_permute2f128_pd(vj0_rel, vk0_rel, 0x31);
    vk0_img = _mm256_permute2f128_pd(vj0_img, vk0_img, 0x31);
    // Second Half
    vt1_rel = _mm256_permute2f128_pd(vj1_rel, vk1_rel, 0x20);
    vt1_img = _mm256_permute2f128_pd(vj1_img, vk1_img, 0x20);
    vk1_rel = _mm256_permute2f128_pd(vj1_rel, vk1_rel, 0x31);
    vk1_img = _mm256_permute2f128_pd(vj1_img, vk1_img, 0x31);
    // GS
    // First Half
    GS_BUTTERFLY(vt0_rel, vt0_img, vk0_rel, vk0_img, vw2_rel, vw2_img, vj0_rel,
                 vj0_img);
    // Second Half
    GS_BUTTERFLY(vt1_rel, vt1_img, vk1_rel, vk1_img, vw2_rel, vw2_img, vj1_rel,
                 vj1_img);
    // Forth Layer (len = 1), four roots per half
    // Transpose
    // First Half
    vj0_rel = _mm256_unpacklo_pd(vt0_rel, vk0_rel);
    vj0_img = _mm256_unpacklo_pd(vt0_img, vk0_img);
    vk0_rel = _mm256_unpackhi_pd(vt0_rel, vk0_rel);
    vk0_img = _mm256_unpackhi_pd(vt0_img, vk0_img);
    // Second Half
    vj1_rel = _mm256_unpacklo_pd(vt1_rel, vk1_rel);
    vj1_img = _mm256_unpacklo_pd(vt1_img, vk1_img);
    vk1_rel = _mm256_unpackhi_pd(vt1_rel, vk1_rel);
    vk1_img = _mm256_unpackhi_pd(vt1_img, vk1_img);
    // First Half
    vw_rel = _mm256_set1_pd(w1_rel);
    vw_img = _mm256_set1_pd(w1_img);
    GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Tranpose Back
    // First Half
    vt0_rel = _mm256_unpacklo_pd(vj0_rel, vk0_rel);
    vt0_img = _mm256_unpacklo_pd(vj0_img, vk0_img);
    vk0_rel = _mm256_unpackhi_pd(vj0_rel, vk0_rel);
    vk0_img = _mm256_unpackhi_pd(vj0_img, vk0_img);
    vj0_rel = _mm256_permute2f128_pd(vt0_rel, vk0_rel, 0x20);
    vj0_img = _mm256_permute2f128_pd(vt0_img, vk0_img, 0x20);
    vk0_rel = _mm256_permute2f128_pd(vt0_rel, vk0_rel, 0x31);
    vk0_img = _mm256_permute2f128_pd(vt0_img, vk0_img, 0x31);
    // Second Half
    vt1_rel = _mm256_unpacklo_pd(vj1_rel, vk1_rel);
    vt1_img = _mm256_unpacklo_pd(vj1_img, vk1_img);
    vk1_rel = _mm256_unpackhi_pd(vj1_rel, vk1_rel);
    vk1_img = _mm256_unpackhi_pd(vj1_img, vk1_img);
    vj1_rel = _mm256_permute2f128_pd(vt1_rel, vk1_rel, 0x20);
    vj1_img = _mm256_permute2f128_pd(vt1_img, vk1_img, 0x20);
    vk1_rel = _mm256_permute2f128_pd(vt1_rel, vk1_rel, 0x31);
    vk1_img = _mm256_permute2f128_pd(vt1_img, vk1_img, 0x31);
    // Store
    // First Half
    _mm256_storeu_pd(&v->rel[j], vj0_rel);
    _mm256_storeu_pd(&v->img[j], vj0_img);
    _mm256_storeu_pd(&v->rel[j + 4], vk0_rel);
    _mm256_storeu_pd(&v->img[j + 4], vk0_img);
    // Second Half
    _mm256_storeu_pd(&v->rel[j + 8], vj1_rel);
    _mm256_storeu_pd(&v->img[j + 8], vj1_img);
    _mm256_storeu_pd(&v->rel[j + 12], vk1_rel);
    _mm256_storeu_pd(&v->img[j + 12], vk1_img);
  }
#else
  for (uint32_t len = m / 2; len > 0; len >>= 1) {
#pragma omp for
    for (uint32_t i = 0; i < m; i += 2 * len) {
      for (uint32_t j = i; j < i + len; j++) {
        double w_rel = w->rel[len + j - i];
        double w_img = -w->img[len + j - i];  // Conjugate
        double j_rel = v->rel[j];
        double j_img = v->img[j];
        double k_rel = v->rel[j + len];
        double k_img = v->img[j + len];
        double t_rel = j_rel - k_rel;
        double t_img = j_img - k_img;
        v->rel[j] = j_rel + k_rel;
        v->img[j] = j_img + k_img;
        v->rel[j + len] = rel_mul(t_rel, t_img, w_rel, w_img);
        v->img[j + len] = img_mul(t_rel, t_img, w_rel, w_img);
      }
    }
  }
#endif
}
