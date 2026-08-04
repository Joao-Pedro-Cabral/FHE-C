
#include "vector.h"

#include <stdlib.h>

#include "complex.h"

void vec_new(vec* v, int n, int int_part) {
  v->rel = malloc(2 * n * sizeof(int64_t));
  v->img = &v->rel[n];
  v->int_part = int_part;
}

void vec_free(vec* v) { free(v->rel); }

void vec_hermitian(vec* v2, const vec* v1, int n) {
  for (int i = 0; i < n / 2; i++) {
    v2->rel[i] = v1->rel[i];
    v2->rel[n - 1 - i] = v1->rel[i];
    v2->img[i] = -v1->img[i];
    v2->img[n - 1 - i] = -v1->img[i];
  }
}

void fft(vec* v, const vec* w, int n) {
  int64_t w_rel, w_img;
  int64_t t_rel, t_img;

  int k = 1;
  for (int len = n / 2; len >= 2; len >>= 1) {
    int j = 0;
    for (int start = 0; start < n; start = j + len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k++;
      for (j = start; j < start + len; j++) {
        // |w| <= 1, then multiplying doesn't increase the integer part
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        // Each add/sub increases the integer part in one bit
        v->rel[j + len] = fixed_sub(v->rel[j], t_rel);
        v->img[j + len] = fixed_sub(v->img[j], t_img);
        v->rel[j] = fixed_add(v->rel[j], t_rel);
        v->img[j] = fixed_add(v->img[j], t_img);
      }
    }
  }
}

void ifft(vec* v, const vec* w, int n) {
  int64_t w_rel, w_img;
  int64_t t_rel, t_img;

  int k = n / 2 - 1;
  for (int len = 2; len <= n / 2; len <<= 1) {
    int j = 0;
    for (int start = 0; start < n; start = j + len) {
      w_rel = w->rel[k];
      w_img = w->img[k];
      k--;
      for (j = start; j < start + len; j++) {
        t_rel = v->rel[j];
        t_img = v->img[j];
        // Each add/sub increases the integer part in one bit
        v->rel[j] = fixed_add(v->rel[j + len], t_rel);
        v->img[j] = fixed_add(v->img[j + len], t_img);
        v->rel[j + len] = fixed_sub(v->rel[j + len], t_rel);
        v->img[j + len] = fixed_sub(v->img[j + len], t_img);
        // |w| <= 1, then multiplying doesn't increase the integer part
        t_rel = rel_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        t_img = img_mul(v->rel[j + len], v->img[j + len], w_rel, w_img);
        v->rel[j + len] = t_rel;
        v->img[j + len] = t_img;
      }
    }
  }
  // The implicit 1/n removes the extra log2(n) bits of the integer part
}
