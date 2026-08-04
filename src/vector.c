
#include "vector.h"

#include <stdlib.h>

#include "complex.h"
#include "utils.h"

void vec_new(vec* v, uint32_t n, uint32_t int_part) {
  v->rel = malloc(2 * n * sizeof(int64_t));
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

// Q1.63 format for w0, 2 <= n <= 2^20
static const int64_t w0_real[20] = {
    INT64_C(565),
    INT64_C(6521908912666391552),
    INT64_C(8521284645587064832),
    INT64_C(9046147529429925888),
    INT64_C(9178958979493304320),
    INT64_C(9212262081300396032),
    INT64_C(9220594129639322624),
    INT64_C(9222677533903514624),
    INT64_C(9223198409482717184),
    INT64_C(9223328629909620736),
    INT64_C(9223361185112102912),
    INT64_C(9223369323918708736),
    INT64_C(9223371358620734464),
    INT64_C(9223371867296264192),
    INT64_C(9223371994465147904),
    INT64_C(9223372026257369088),
    INT64_C(9223372034205423616),
    INT64_C(9223372036192438272),
    INT64_C(9223372036689190912),
    INT64_C(9223372036813379584),
};

static const int64_t w0_imag[20] = {
    INT64_C(9223372036854775807), INT64_C(6521908912666390528),
    INT64_C(3529631669043774976), INT64_C(1799390620744550912),
    INT64_C(904048551248139776),  INT64_C(452569415305004288),
    INT64_C(226352880908292032),  INT64_C(113184963073701952),
    INT64_C(56593546894377864),   INT64_C(28296906617819880),
    INT64_C(14148469955268186),   INT64_C(7074237058429792),
    INT64_C(3537118789314387),    INT64_C(1768559427169631),
    INT64_C(884279717648870),     INT64_C(442139859332442),
    INT64_C(221069929729722),     INT64_C(110534964872798),
    INT64_C(55267482437391),      INT64_C(27633741218820),
};

void gen_w(vec* w, uint32_t n) {
  vec tmp;
  vec_new(&tmp, n, 1);
  uint32_t index = floor_log2(n);
  int64_t w0_rel = w0_real[index];
  int64_t w0_img = w0_imag[index];
  tmp.rel[0] = (1UL << 63) - 1;
  tmp.img[0] = 0;
  tmp.rel[1] = w0_rel;
  tmp.img[1] = w0_img;
  int64_t w02_rel = rel_mul(w0_rel, w0_img, w0_rel, w0_img);
  int64_t w02_img = img_mul(w0_rel, w0_img, w0_rel, w0_img);
  // Compute in normal order (Q1.63 format)
  for (uint32_t i = 2; i < n; i++) {
    tmp.rel[i] = rel_mul(tmp.rel[i - 1], tmp.img[i - 1], w02_rel, w02_img);
    tmp.img[i] = img_mul(tmp.rel[i - 1], tmp.img[i - 1], w02_rel, w02_img);
  }
  // Then change to bit-reverse order
  for (uint32_t i = 1; i < n; i++) {
    uint32_t j = bit_reverse32(i, index);
    w->rel[j] = tmp.rel[i];
    w->img[j] = tmp.img[i];
  }
  vec_free(&tmp);
}

void fft(vec* v, const vec* w, uint32_t n) {
  int64_t w_rel, w_img;
  int64_t t_rel, t_img;

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
        v->rel[j + len] = fixed_sub(v->rel[j], t_rel);
        v->img[j + len] = fixed_sub(v->img[j], t_img);
        v->rel[j] = fixed_add(v->rel[j], t_rel);
        v->img[j] = fixed_add(v->img[j], t_img);
      }
    }
  }
}

void ifft(vec* v, const vec* w, uint32_t n) {
  int64_t w_rel, w_img;
  int64_t t_rel, t_img;

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
