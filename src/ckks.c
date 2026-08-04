
#include "ckks.h"

#include <string.h>

void ckks_keygen(vec* w, int n) {
  gen_w(w, n);
}

void ckks_encode(poly* scaled_m, const vec* v, const vec* w, uint64_t scale,
                 int n) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, n, v->int_part);
  vec_hermitian(&v2, v, n);
  // Second step: IFFT
  ifft(&v2, w, n);
  // Third step: Rescale
  // They're equal, because X^n + 1 only has real coeffs
  poly* m = v2.rel;
  poly_rescale(scaled_m, m, scale, v2.int_part, n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly* scaled_m, const vec* w, uint64_t scale, int n) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, n, v->int_part);
  // They're equal, because X^n + 1 only has real coeffs
  poly* m = v2.rel;
  poly_rescale(m, scaled_m, scale, v2.int_part, n);
  memset(v2.img, 0, n);
  // Second step: FFT
  fft(&v2, w, n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, n);
  // Free
  vec_free(&v2);
}
