
#include "rlwe.h"

#include <stdlib.h>
#include <string.h>

void rlwe_new(rlwe* a, uint32_t n, uint32_t l) {
  *a = malloc(sizeof(uint64_t*) * l);
  for (uint32_t i = 0; i < l; i++) {
    poly_new(&(*a)[i], n);
  }
}

void rlwe_free(rlwe* a, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_free(&(*a)[i]);
  }
  free(*a);
}

void rlwe_gen_secret(rlwe s, const set eta, const uint64_t* q, uint32_t n,
                     uint32_t l) {
  poly s_hat;
  poly_new(&s_hat, n);
  poly_ternary_sampler(s_hat, n);
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s[i], s_hat, sizeof(int64_t) * n);
    poly_ntt(s[i], eta[i], n);
  }
  poly_free(&s_hat);
}
