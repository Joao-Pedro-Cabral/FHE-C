
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"
#include "random.h"
#include "xof.h"

void test_cbd(poly_t e, uint32_t n) {
  uint8_t *buf = malloc(6 * n * sizeof(uint8_t));
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  xof(buf, 6 * n, seed, SEED_BYTES);
  printf("0x");
  for (uint32_t i = 0; i < 6 * n; i++) {
    printf("%02x", buf[i]);
  }
  printf("\n");
  uint8_t *b0_8 = buf + 0 * n, *b1_8 = buf + 1 * n, *b2_8 = buf + 2 * n;
  uint8_t *c0_8 = buf + 3 * n, *c1_8 = buf + 4 * n, *c2_8 = buf + 5 * n;
  for (uint32_t i = 0; i < n; i++) {
    int b = __builtin_popcount(b0_8[i]) + __builtin_popcount(b1_8[i]) +
            __builtin_popcount(b2_8[i]);
    int c = __builtin_popcount(c0_8[i]) + __builtin_popcount(c1_8[i]) +
            __builtin_popcount(c2_8[i]);
    e[i] = b - c;
  }
  free(buf);
}

int main() {
  poly_t e1, e2;
  uint32_t n = 16;
  poly_new(&e1, n);
  poly_new(&e2, n);

  poly_cbd24_sampler(e1, n);
  test_cbd(e2, n);

  for (uint32_t i = 0; i < n; i++) {
    if (e1[i] != e2[i]) {
      printf("e1[%d]: %ld, e2[%d]: %ld\n", i, e1[i], i, e2[i]);
    }
  }

  poly_free(&e1);
  poly_free(&e2);

  return 0;
}
