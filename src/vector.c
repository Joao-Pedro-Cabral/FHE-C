
#include "vector.h"

#include <math.h>
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
    v2->img[i] = -v1->img[i];
    v2->img[n - 1 - i] = -v1->img[i];
  }
}

void vec_input(vec* v2, const vec* v1, uint32_t n) {
  for (uint32_t i = 0; i < n / 2; i++) {
    v2->rel[i] = v1->rel[i];
    v2->img[i] = 0;
  }
}

void gen_w(vec* w, uint32_t n) {
  vec tmp;
  vec_new(&tmp, n, 1);
  uint32_t log2_n = floor_log2(n);
  // 1st quadrant
  for (uint32_t i = 1; i < n / 2; i += 2) {
    w->rel[bit_reverse32((i - 1) / 2, log2_n)] = cos(i * M_PI / n);
    w->img[bit_reverse32((i - 1) / 2, log2_n)] = sin(i * M_PI / n);
  }
  // 2nd quadrant
  for (uint32_t i = 1; i < n / 2; i += 2) {
    w->rel[bit_reverse32((i + n / 2 - 1) / 2, log2_n)] =
        -w->rel[bit_reverse32((i - 1) / 2, log2_n)];
    w->img[bit_reverse32((i + n / 2 - 1) / 2, log2_n)] =
        w->img[bit_reverse32((i - 1) / 2, log2_n)];
  }
  // 3rd quadrant
  for (uint32_t i = 1; i < n / 2; i += 2) {
    w->rel[bit_reverse32((i + n - 1) / 2, log2_n)] =
        -w->rel[bit_reverse32((i - 1) / 2, log2_n)];
    w->img[bit_reverse32((i + n - 1) / 2, log2_n)] =
        -w->img[bit_reverse32((i - 1) / 2, log2_n)];
  }
  // 4th quadrant
  for (uint32_t i = 1; i < n / 2; i += 2) {
    w->rel[bit_reverse32((i + 3 * n / 2 - 1) / 2, log2_n)] =
        w->rel[bit_reverse32((i - 1) / 2, log2_n)];
    w->img[bit_reverse32((i + 3 * n / 2 - 1) / 2, log2_n)] =
        -w->img[bit_reverse32((i - 1) / 2, log2_n)];
  }
  vec_free(&tmp);
}

void fft(vec* v, const vec* w, uint32_t n) {
  double w_rel, w_img;
  double t_rel, t_img;

  uint32_t k = 1;
  for (uint32_t len = n / 2; len >= 2; len >>= 1) {
    uint32_t j = 0;
    for (uint32_t start = 0; start < n; start = j + len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k++;
      for (j = start; j < start + len; j++) {
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

  uint32_t k = n / 2 - 1;
  for (uint32_t len = 2; len <= n / 2; len <<= 1) {
    uint32_t j = 0;
    for (uint32_t start = 0; start < n; start = j + len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k--;
      for (j = start; j < start + len; j++) {
        t_rel = v->rel[j];
        t_img = v->img[j];
        // Each add/sub increases the integer part in one bit
        v->rel[j] = v->rel[j + len] + t_rel;
        v->img[j] = v->img[j + len] + t_img;
        v->rel[j + len] = v->rel[j + len] - t_rel;
        v->img[j + len] = v->img[j + len] - t_img;
        // |w| <= 1, then multiplying doesn't increase the integer part
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        v->rel[j + len] = t_rel;
        v->img[j + len] = t_img;
      }
    }
  }
  // The implicit 1/n removes the extra log2(n) bits of the integer part
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = v->rel[i] / n;
    v->img[i] = v->img[i] / n;
  }
}
