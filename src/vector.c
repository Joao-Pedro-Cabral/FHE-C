
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

void ifft(vec* v, const vec* w, uint32_t n) {
  double w_rel, w_img;
  double t_rel, t_img;

  for (uint32_t len = 1; len < n; len <<= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      // We only have the upper half roots, so we conjugate here
      w_rel = w->rel[k];
      w_img = -w->img[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t_rel = v->rel[j];
        t_img = v->img[j];
        v->rel[j] = v->rel[j + len] + t_rel;
        v->img[j] = v->img[j + len] + t_img;
        v->rel[j + len] = t_rel - v->rel[j + len];
        v->img[j + len] = t_img - v->img[j + len];
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        v->rel[j + len] = t_rel;
        v->img[j + len] = t_img;
      }
    }
  }
  // Last layer (len = 1), everyone uses the same root
  // __m256d vw_rel = _mm256_set1_pd(w->rel[1]);
  // __m256d vw_img = _mm256_set1_pd(-w->img[1]);
  // for (uint32_t j = 0; j < n; j += 8) {
  //   // Load
  //   __m256d vj_rel = _mm256_loadu_pd(&v->rel[j]);
  //   __m256d vj_img = _mm256_loadu_pd(&v->img[j]);
  //   __m256d vk_rel = _mm256_loadu_pd(&v->rel[j + 4]);
  //   __m256d vk_img = _mm256_loadu_pd(&v->img[j + 4]);
  //   // Tranpose
  //   __m256d vt_rel = _mm256_unpacklo_pd(vj_rel, vk_rel);
  //   __m256d vt_img = _mm256_unpacklo_pd(vj_img, vk_img);
  //   vk_rel = _mm256_unpackhi_pd(vj_rel, vk_rel);
  //   vk_img = _mm256_unpackhi_pd(vj_img, vk_img);
  //   // GS
  //   vj_rel = _mm256_add_pd(vt_rel, vk_rel);
  //   vj_img = _mm256_add_pd(vt_img, vk_img);
  //   vk_rel = _mm256_sub_pd(vt_rel, vk_rel);
  //   vk_img = _mm256_sub_pd(vt_img, vk_img);
  //   vt_rel = rel_vmul(vk_rel, vk_img, vw_rel, vw_img);
  //   vt_img = img_vmul(vk_rel, vk_img, vw_rel, vw_img);
  //   // Transpose back
  //   vk_rel = _mm256_unpackhi_pd(vj_rel, vt_rel);
  //   vk_img = _mm256_unpackhi_pd(vj_img, vt_img);
  //   vj_rel = _mm256_unpacklo_pd(vj_rel, vt_rel);
  //   vj_img = _mm256_unpacklo_pd(vj_img, vt_img);
  //   // Store
  //   _mm256_storeu_pd(&v->rel[j], vj_rel);
  //   _mm256_storeu_pd(&v->img[j], vj_img);
  //   _mm256_storeu_pd(&v->rel[j + 4], vk_rel);
  //   _mm256_storeu_pd(&v->img[j + 4], vk_img);
  // }
  double inv_n = 1.0 / n;
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = v->rel[i] * inv_n;
    v->img[i] = v->img[i] * inv_n;
  }
}
