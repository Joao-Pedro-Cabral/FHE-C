
#include "vector.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "complex.h"
#include "utils.h"

void vec_new(vec* v, uint32_t n, uint32_t int_part) {
  v->rel = malloc(2 * n * sizeof(double));
  v->img = &v->rel[n];
  v->int_part = int_part;
}

void vec_free(vec* v) { free(v->rel); }

void vec_hermitian(vec* v2, const vec* v1, uint32_t n) {
  for (uint32_t i = 0; i < n / 2; i++) {
    v2->rel[i] = v1->rel[i];
    v2->rel[n - 1 - i] = v1->rel[i];
    v2->img[i] = v1->img[i];
    v2->img[n - 1 - i] = -v1->img[i];
  }
}

void vec_input(vec* v2, const vec* v1, uint32_t n) {
  for (uint32_t i = 0; i < n / 2; i++) {
    v2->rel[i] = v1->rel[i];
    v2->img[i] = 0;
  }
}

void vec_gen_w(vec* w, uint32_t n) {
  uint32_t log2_2n = floor_log2(2 * n);
  // 1st quadrant
  w->rel[bit_reverse32(0, log2_2n)] = 1.0;
  w->img[bit_reverse32(0, log2_2n)] = 0.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(i, log2_2n)] = cos(i * M_PI / n);
    w->img[bit_reverse32(i, log2_2n)] = sin(i * M_PI / n);
  }
  // 2nd quadrant
  w->rel[bit_reverse32(n / 2, log2_2n)] = 0.0;
  w->img[bit_reverse32(n / 2, log2_2n)] = 1.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(i + n / 2, log2_2n)] =
        -w->rel[bit_reverse32(i, log2_2n)];
    w->img[bit_reverse32(i + n / 2, log2_2n)] =
        w->img[bit_reverse32(i, log2_2n)];
  }
  // 3rd quadrant
  w->rel[bit_reverse32(n, log2_2n)] = -1.0;
  w->img[bit_reverse32(n, log2_2n)] = 0.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(i + n, log2_2n)] = -w->rel[bit_reverse32(i, log2_2n)];
    w->img[bit_reverse32(i + n, log2_2n)] = -w->img[bit_reverse32(i, log2_2n)];
  }
  // 4th quadrant
  w->rel[bit_reverse32(3 * n / 2, log2_2n)] = 0.0;
  w->img[bit_reverse32(3 * n / 2, log2_2n)] = -1.0;
  for (uint32_t i = 1; i < n / 2; i++) {
    w->rel[bit_reverse32(i + 3 * n / 2, log2_2n)] =
        w->rel[bit_reverse32(i, log2_2n)];
    w->img[bit_reverse32(i + 3 * n / 2, log2_2n)] =
        -w->img[bit_reverse32(i, log2_2n)];
  }
}

void fft(vec* v, const vec* w, uint32_t n) {
  double w_rel, w_img;
  double t_rel, t_img;

  uint32_t k = 1;
  for (uint32_t len = n / 2; len > 0; len >>= 1) {
    for (uint32_t start = 0; start < n; start += 2 * len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k++;
      for (uint32_t j = start; j < start + len; j++) {
        // |w| <= 1, then multiplying doesn't increase the integer part
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        // Each add/sub increases the integer part in one bit
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

  uint32_t k = floor_log2(n);
  for (uint32_t l = 0; l < k; l++) {
    for (uint32_t i = 0; i < (1UL << (k - 1 - l)); i++) {
      w_rel = w->rel[(1UL << (k - 1 - l)) + i];
      w_img = -w->img[(1UL << (k - 1 - l)) + i];
      uint32_t len = (1UL << l);
      for (uint32_t j = i * (1UL << (l + 1)); j < i * (1UL << (l + 1)) + len;
           j++) {
        t_rel = v->rel[j];
        t_img = v->img[j];
        // Each add/sub increases the integer part in one bit
        v->rel[j] = v->rel[j + len] + t_rel;
        v->img[j] = v->img[j + len] + t_img;
        v->rel[j + len] = t_rel - v->rel[j + len];
        v->img[j + len] = t_img - v->img[j + len];
        // |w| <= 1, then multiplying doesn't increase the integer part
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        v->rel[j + len] = t_rel;
        v->img[j + len] = t_img;
      }
    }
  }
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = v->rel[i] / n;
    v->img[i] = v->img[i] / n;
  }
}
