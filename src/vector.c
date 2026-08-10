
#include "vector.h"

#include <immintrin.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "complex.h"
#include "utils.h"

void vec_new(vec* v, uint32_t n) {
  v->rel = malloc(2 * n * sizeof(double));
  v->img = &v->rel[n];
}

void vec_free(vec* v) { free(v->rel); }

// TODO: Vectorize only in forward mode
void vec_hermitian(vec* v2, const vec* v1, uint32_t n) {
  uint32_t log_n = floor_log2(n);
  // NO -> BO for IFFT
  for (uint32_t i = 0; i < n / 2; i++) {
    v2->rel[bit_reverse32(i, log_n)] = v1->rel[i];
    v2->rel[bit_reverse32(n - 1 - i, log_n)] = v1->rel[i];
    v2->img[bit_reverse32(i, log_n)] = v1->img[i];
    v2->img[bit_reverse32(n - 1 - i, log_n)] = -v1->img[i];
  }
}

void vec_input(vec* v2, const vec* v1, uint32_t n) {
  uint32_t log_n = floor_log2(n);
  // BO -> NO for FFT
  for (uint32_t i = 0; i < n / 2; i++) {
    v2->rel[i] = v1->rel[bit_reverse32(i, log_n)];
    v2->img[i] = v1->img[bit_reverse32(i, log_n)];
  }
}

void vec_gen_roots(vec* w, uint32_t n) {
  uint32_t log2_n = floor_log2(n);
  // 1st quadrant
  w->rel[bit_reverse32(0, log2_n)] = 1.0;
  w->img[bit_reverse32(0, log2_n)] = 0.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(i, log2_n)] = cos(i * M_PI / n);
    w->img[bit_reverse32(i, log2_n)] = sin(i * M_PI / n);
  }
  // 2nd quadrant
  w->rel[bit_reverse32(n / 2, log2_n)] = 0.0;
  w->img[bit_reverse32(n / 2, log2_n)] = 1.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(n - i, log2_n)] = -w->rel[bit_reverse32(i, log2_n)];
    w->img[bit_reverse32(n - i, log2_n)] = w->img[bit_reverse32(i, log2_n)];
  }
}

void fft(vec* v, const vec* w, uint32_t n) {
  double w_rel, w_img;
  double t_rel, t_img;

  for (uint32_t len = n / 2; len > 0; len >>= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        v->rel[j + len] = v->rel[j] - t_rel;
        v->img[j + len] = v->img[j] - t_img;
        v->rel[j] = v->rel[j] + t_rel;
        v->img[j] = v->img[j] + t_img;
      }
    }
  }
}

#define GS_BUTTERFLY(vj_rel, vj_img, vk_rel, vk_img, vw_rel, vw_img, vt_rel, \
                     vt_img)                                                 \
  do {                                                                       \
    vt_rel = _mm256_sub_pd((vj_rel), (vk_rel));                              \
    vt_img = _mm256_sub_pd((vj_img), (vk_img));                              \
    (vj_rel) = _mm256_add_pd((vj_rel), (vk_rel));                            \
    (vj_img) = _mm256_add_pd((vj_img), (vk_img));                            \
    (vk_rel) = rel_vmul(vt_rel, vt_img, (vw_rel), (vw_img));                 \
    (vk_img) = img_vmul(vt_rel, vt_img, (vw_rel), (vw_img));                 \
  } while (0)

void ifft(vec* v, const vec* w, uint32_t n) {
  uint32_t k1 = n / 2;
  uint32_t k2 = n / 4;
  uint32_t k3 = n / 8;
  uint32_t k4 = n / 16;
  __m256d sign = _mm256_set1_pd(-0.0f);
  // Merge 4 layers using all the 16 AVX2 registers
  for (uint32_t j = 0; j < n; j += 16) {
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
    __m256d vw_rel = _mm256_loadu_pd(&w->rel[k1]);
    __m256d vw_img = _mm256_loadu_pd(&w->img[k1]);
    vw_img = _mm256_xor_pd(vw_img, sign);
    k1 += 4;
    GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    vw_rel = _mm256_loadu_pd(&w->rel[k1]);
    vw_img = _mm256_loadu_pd(&w->img[k1]);
    vw_img = _mm256_xor_pd(vw_img, sign);
    k1 += 4;
    GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
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
    vw_rel =
        _mm256_setr_pd(w->rel[k2], w->rel[k2], w->rel[k2 + 1], w->rel[k2 + 1]);
    vw_img = _mm256_setr_pd(-w->img[k2], -w->img[k2], -w->img[k2 + 1],
                            -w->img[k2 + 1]);
    k2 += 2;
    GS_BUTTERFLY(vt0_rel, vt0_img, vk0_rel, vk0_img, vw_rel, vw_img, vj0_rel,
                 vj0_img);
    // Second Half
    vw_rel =
        _mm256_setr_pd(w->rel[k2], w->rel[k2], w->rel[k2 + 1], w->rel[k2 + 1]);
    vw_img = _mm256_setr_pd(-w->img[k2], -w->img[k2], -w->img[k2 + 1],
                            -w->img[k2 + 1]);
    k2 += 2;
    GS_BUTTERFLY(vt1_rel, vt1_img, vk1_rel, vk1_img, vw_rel, vw_img, vj1_rel,
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
    vw_rel = _mm256_set1_pd(w->rel[k3]);
    vw_img = _mm256_set1_pd(-w->img[k3]);
    k3 += 1;
    GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    // Second Half
    vw_rel = _mm256_set1_pd(w->rel[k3]);
    vw_img = _mm256_set1_pd(-w->img[k3]);
    k3 += 1;
    GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
                 vt1_img);
    // Forth layer (len = 8), one root
    vw_rel = _mm256_set1_pd(w->rel[k4]);
    vw_img = _mm256_set1_pd(-w->img[k4]);
    k4 += 1;
    GS_BUTTERFLY(vj0_rel, vj0_img, vj1_rel, vj1_img, vw_rel, vw_img, vt0_rel,
                 vt0_img);
    GS_BUTTERFLY(vk0_rel, vk0_img, vk1_rel, vk1_img, vw_rel, vw_img, vt1_rel,
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
  uint32_t log2_n = floor_log2(n);
  uint32_t n4 = ((log2_n & 0x1) == 0) ? n : (n / 2);
  for (uint32_t len = 16; len < n4; len <<= 2) {
    uint32_t klen = n / (2 * len);
    uint32_t k2len = n / (4 * len);
    for (uint32_t start = 0; start < n; start += 4 * len) {
      // We only have the upper half roots, so we conjugate here
      // First Layer: two roots
      __m256d vw10_rel = _mm256_set1_pd(w->rel[klen]);
      __m256d vw10_img = _mm256_set1_pd(-w->img[klen]);
      __m256d vw11_rel = _mm256_set1_pd(w->rel[klen + 1]);
      __m256d vw11_img = _mm256_set1_pd(-w->img[klen + 1]);
      klen += 2;
      // Second Layer: one root
      __m256d vw2_rel = _mm256_set1_pd(w->rel[k2len]);
      __m256d vw2_img = _mm256_set1_pd(-w->img[k2len]);
      k2len += 1;
      for (uint32_t j = start; j < start + len; j += 4) {
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
        // First Layer (len)
        __m256d vt_rel, vt_img;
        GS_BUTTERFLY(vj0_rel, vj0_img, vk0_rel, vk0_img, vw10_rel, vw10_img,
                     vt_rel, vt_img);
        GS_BUTTERFLY(vj1_rel, vj1_img, vk1_rel, vk1_img, vw11_rel, vw11_img,
                     vt_rel, vt_img);
        // Second layer (2*len)
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
        _mm256_storeu_pd(&v->rel[j + 2 * len], vj1_rel);
        _mm256_storeu_pd(&v->img[j + 2 * len], vj1_img);
        _mm256_storeu_pd(&v->rel[j + 3 * len], vk1_rel);
        _mm256_storeu_pd(&v->img[j + 3 * len], vk1_img);
      }
    }
  }

  // Final layer for odd log2_n
  if (n4 < n) {
    // Only one root for everyone
    __m256d vw_rel = _mm256_set1_pd(w->rel[1]);
    __m256d vw_img = _mm256_set1_pd(-w->img[1]);
    // Unroll two vector iterations
    for (uint32_t j = 0; j < n / 2; j += 8) {
      // Load
      // First Half
      __m256d vj0_rel = _mm256_loadu_pd(&v->rel[j]);
      __m256d vj0_img = _mm256_loadu_pd(&v->img[j]);
      __m256d vk0_rel = _mm256_loadu_pd(&v->rel[j + n / 2]);
      __m256d vk0_img = _mm256_loadu_pd(&v->img[j + n / 2]);
      // Second Half
      __m256d vj1_rel = _mm256_loadu_pd(&v->rel[j + 4]);
      __m256d vj1_img = _mm256_loadu_pd(&v->img[j + 4]);
      __m256d vk1_rel = _mm256_loadu_pd(&v->rel[j + 4 + n / 2]);
      __m256d vk1_img = _mm256_loadu_pd(&v->img[j + 4 + n / 2]);
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
      _mm256_storeu_pd(&v->rel[j + n / 2], vk0_rel);
      _mm256_storeu_pd(&v->img[j + n / 2], vk0_img);
      // Second Half
      _mm256_storeu_pd(&v->rel[j + 4], vj1_rel);
      _mm256_storeu_pd(&v->img[j + 4], vj1_img);
      _mm256_storeu_pd(&v->rel[j + 4 + n / 2], vk1_rel);
      _mm256_storeu_pd(&v->img[j + 4 + n / 2], vk1_img);
    }
  }

  double inv_n = 1.0 / n;
  __m256d vinvn = _mm256_set1_pd(inv_n);
  for (uint32_t i = 0; i < n; i += 4) {
    __m256d vj_rel = _mm256_loadu_pd(&v->rel[i]);
    __m256d vj_img = _mm256_loadu_pd(&v->img[i]);
    vj_rel = _mm256_mul_pd(vj_rel, vinvn);
    vj_img = _mm256_mul_pd(vj_img, vinvn);
    _mm256_storeu_pd(&v->rel[i], vj_rel);
    _mm256_storeu_pd(&v->img[i], vj_img);
  }
}
