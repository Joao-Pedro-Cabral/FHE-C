
#include "ckks.h"

#include <string.h>

void ckks_keygen(vec* w, uint32_t n) { gen_w(w, n); }

void ckks_encode(poly m, const vec* v, const vec* w, int64_t scale,
                 uint32_t n) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, n, v->int_part);
  vec_hermitian(&v2, v, n);
  // Second step: IFFT
  ifft(&v2, w, n);
  // Third step: Rescale
  vector_to_poly(m, &v2, scale, n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly m, const vec* w, int64_t scale,
                 uint32_t n) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, n, v->int_part);
  poly_to_vector(&v2, m, scale, n);
  // Second step: FFT
  fft(&v2, w, n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, n);
  // Free
  vec_free(&v2);
}
