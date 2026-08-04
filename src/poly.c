
#include "poly.h"

#include "fixed.h"

void poly_new(poly* m, int n) {
  m = malloc(n * sizeof(int64_t));
}

void poly_free(poly* m) { free(m); }

void poly_rescale(poly* m2, const poly* m, uint64_t scale, int int_part,
                  int n) {
  // Since >> 63 = >> (int_part = frac_part)
  // Then, we multiply scale by (1UL << int_part)
  // For completely eliminate the fractional part with a multiplication
  scale = scale * (1UL << int_part);
  for (int i = 0; i < n; i++) {
    m2[i] = fixed_mul(m[i], scale);
  }
}
