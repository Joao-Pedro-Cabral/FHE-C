
#include "vector.h"

#include <immintrin.h>
#include <math.h>
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

// m = n/2 (half cyclotomic size)
void fft(vec_t* v, const vec_t* w, uint32_t m) {
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
        __m256d vw21_rel = _mm256_loadu_pd(&w->rel[k2len + j + 4]);
        __m256d vw21_img = _mm256_loadu_pd(&w->img[k2len + j + 4]);
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
      // Four roots
      __m256d vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j]);
      __m256d vw_img = _mm256_loadu_pd(&w->img[m / 2 + j]);
      // CT
      __m256d vt_rel, vt_img;
      CT_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
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
}

#define GS_BUTTERFLY(vj_rel, vj_img, vk_rel, vk_img, vw_rel, vw_img, vt_rel, \
                     vt_img)                                                 \
  do {                                                                       \
    (vt_rel) = _mm256_sub_pd((vj_rel), (vk_rel));                            \
    (vt_img) = _mm256_sub_pd((vj_img), (vk_img));                            \
    (vj_rel) = _mm256_add_pd((vj_rel), (vk_rel));                            \
    (vj_img) = _mm256_add_pd((vj_img), (vk_img));                            \
    (vk_rel) = rel_vmul((vt_rel), (vt_img), (vw_rel), (vw_img));             \
    (vk_img) = img_vmul((vt_rel), (vt_img), (vw_rel), (vw_img));             \
  } while (0)

// m = n/2 (half cyclotomic size)
void ifft(vec_t* v, const vec_t* w, uint32_t m) {
  uint32_t log2_m = floor_log2(m);
  uint32_t m2 = ((log2_m & 0x1) == 0) ? m : (m / 2);
  __m256d sign = _mm256_set1_pd(-0.0f);
  // First layer for odd log2_m
  if (m2 < m) {
    // Unroll two vector iterations
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
      // Four roots
      __m256d vw_rel = _mm256_loadu_pd(&w->rel[m / 2 + j]);
      __m256d vw_img = _mm256_xor_pd(_mm256_loadu_pd(&w->img[m / 2 + j]), sign);
      // GS
      __m256d vt_rel, vt_img;
      GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt_rel,
                   vt_img);
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
        __m256d vw11_rel = _mm256_loadu_pd(&w->rel[klen + j + 4]);
        __m256d vw11_img =
            _mm256_xor_pd(_mm256_loadu_pd(&w->img[klen + j + 4]), sign);
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
}
