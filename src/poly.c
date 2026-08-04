
#include "poly.h"

#include <stdlib.h>

#include "fixed.h"

void poly_new(poly* m, uint32_t n) { *m = malloc(n * sizeof(int64_t)); }

void poly_free(poly* m) { free(*m); }

void vector_to_poly(poly m, const vec* v, int64_t scale, uint32_t n) {
  // Since >> 63 = >> (int_part + frac_part)
  // Then, we multiply scale by (1UL << int_part)
  // For completely eliminate the fractional part with a multiplication
  scale = scale * (INT64_C(1) << v->int_part);
  for (uint32_t i = 0; i < n; i++) {
    m[i] = fixed_mul(v->rel[i], scale);
  }
}

void poly_to_vector(vec* v, const poly m, int64_t scale, uint32_t n) {
  // To increase the precision, we first compute (1 << int_part)/scale
  // Instead of /scale and *(1 << int_part)
  scale = ((INT64_C(1) << (63 - v->int_part)) + scale / 2) / scale;
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = m[i] * scale;
    v->img[i] = 0;
  }
}
