
#include "vector.h"

#include <math.h>
#include <stdio.h>
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
    v2->img[i] = v1->img[i];
  }
}

void vec_gen_w(vec* w, uint32_t n) {
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
    w->rel[bit_reverse32(i + n / 2, log2_n)] =
        -w->rel[bit_reverse32(i, log2_n)];
    w->img[bit_reverse32(i + n / 2, log2_n)] = w->img[bit_reverse32(i, log2_n)];
  }
}

void fft(vec* v, const vec* w, uint32_t n) {
  double w_rel, w_img;
  double t_rel, t_img;

  for (uint32_t len = 1; len < n; len <<= 1) {
    uint32_t k = 0;
    uint32_t ratio = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      for (uint32_t j = start; j < start + len; ++j) {
        uint32_t idx = len + k / ratio;
        w_rel = w->rel[idx];
        w_img = w->img[idx];
        k++;
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

  for (uint32_t len = n / 2; len > 0; len >>= 1) {
    uint32_t k = 0;
    uint32_t ratio = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      // We only have the upper half roots, so we conjugate here
      for (uint32_t j = start; j < start + len; ++j) {
        uint32_t idx = len + k / ratio;
        w_rel = w->rel[idx];
        w_img = -w->img[idx];
        k++;
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
  double inv_n = 1.0 / n;
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = v->rel[i] * inv_n;
    v->img[i] = v->img[i] * inv_n;
  }
}
