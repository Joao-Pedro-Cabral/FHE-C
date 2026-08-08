
#include "rns.h"

#include <stdlib.h>
#include <string.h>

void rns_new(rns* a, uint32_t n, uint32_t l) {
  *a = malloc(sizeof(uint64_t*) * l);
  for (uint32_t i = 0; i < l; i++) {
    poly_new(&(*a)[i], n);
  }
}

void rns_free(rns* a, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_free(&(*a)[i]);
  }
  free(*a);
}

void rns_gen_secret_key(rns s, const rns roots, const uint64_t* q, uint32_t n,
                        uint32_t l) {
  poly s_hat;
  poly_new(&s_hat, n);
  poly_ternary_sampler(s_hat, n);
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s[i], s_hat, sizeof(int64_t) * n);
    poly_ntt(s[i], roots[i], q[i], n);
  }
  poly_free(&s_hat);
}

void rns_gen_roots(rns roots, const uint64_t* q, uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_roots(roots[i], n, (int64_t)q[i]);
  }
}

void rns_gen_inverse_roots(rns roots, const uint64_t* q, uint32_t n,
                           uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_inverse_roots(roots[i], n, (int64_t)q[i]);
  }
}
