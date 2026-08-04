
#include "poly.h"

#include <stdlib.h>

#include "fixed.h"

void poly_new(poly* m, uint32_t n) {
  *m = malloc(n * sizeof(int64_t));
}

void poly_free(poly *m) { free(*m); }

void poly_rescale(poly m2, const poly m, int64_t scale, uint32_t int_part,
                  uint32_t n) {
  // Since >> 63 = >> (int_part = frac_part)
  // Then, we multiply scale by (1UL << int_part)
  // For completely eliminate the fractional part with a multiplication
  scale = scale * (INT64_C(1) << int_part);
  for (uint32_t i = 0; i < n; i++) {
    m2[i] = fixed_mul(m[i], scale);
  }
}
