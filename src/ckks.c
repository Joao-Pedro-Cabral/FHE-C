
#include "ckks.h"

#include <string.h>

void ckks_encode(poly* scaled_m, const vec* v, const vec* w, uint64_t scale,
                 int n) {
  // First step: Convert v into v2 = (v, v)
  vec v2;
  vec_new(&v2, n, v->int_part);
  vec_hermitian(&v2, &v);
  // Second step: IFFT
  ifft(&v2, w, n);
  // Third step: Rescale
  // They're equal, because X^n + 1 only has real coeffs
  poly* m = v2.rel;
  poly_rescale(scaled_m, m, scale, v2.int_part, n);
}
