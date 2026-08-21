
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

// m = n/2 (half cyclotomic size)
void fft(vec_t* v, const vec_t* w, uint32_t m) {
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
}

// m = n/2 (half cyclotomic size)
void ifft(vec_t* v, const vec_t* w, uint32_t m) {
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
}
