
#include "set.h"

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

void set_new(set* a, uint32_t n, uint32_t l) {
  *a = malloc(sizeof(roots) * l);
  for (uint32_t i = 0; i < l; i++) {
    (*a)[i] = malloc(sizeof(int64_t) * n);
  }
}

void set_free(set* a, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    free((*a)[i]);
  }
  free(*a);
}

void gen_eta(roots eta, uint32_t n, int64_t q) {
  // Find one primitive root (n > 2)
  uint32_t log2_n = floor_log2(n);
  int64_t eta0;  // our primitive root
  for (eta0 = 2; eta0 < q - 2; eta0++) {
    int64_t x = eta0;
    for (uint32_t j = 0; j < log2_n; j++) {
      x = (x * x) % q;
      if (x > q / 2) x -= q;
    }
    if (x != -1) continue;  // x^(n) = -1
    x = (x * x) % q;
    if (x == 1) break;  // x^(2n) = 1
  }
  eta[0] = 1;
  for (uint32_t i = 0; i < n; i++) {
    eta[bit_reverse32(i, log2_n)] =
        (eta[bit_reverse32(i - 1, log2_n)] * eta0) % q;
  }
}

void set_gen_eta(set eta, const uint64_t* q, uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    gen_eta(eta[i], n, (int64_t)q[i]);
  }
}
