
#include "rns.h"

#include <stdlib.h>
#include <string.h>

#include "mod.h"

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

void rns_gen_secret_key(rns s, const rns roots, const uint64_t* q,
                        const int64_t* qR, const int64_t* qinv, uint32_t n,
                        uint32_t l) {
  poly s_hat;
  poly_new(&s_hat, n);
  poly_ternary_sampler(s_hat, n);
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s[i], s_hat, sizeof(int64_t) * n);
    poly_ntt(s[i], roots[i], q[i], qR[i], qinv[i], n);
  }
  poly_free(&s_hat);
}

void rns_gen_roots(rns roots, const uint64_t* q, const int64_t* qinv,
                   uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_roots(roots[i], n, (int64_t)q[i], qinv[i]);
  }
}

void rns_gen_inverse_roots(rns roots, const uint64_t* q, const int64_t* qinv,
                           uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_inverse_roots(roots[i], n, (int64_t)q[i], qinv[i]);
  }
}

void rns_gen_q_new(int64_t** qn, int64_t** qR, int64_t** qinv, uint32_t l) {
  *qn = malloc(sizeof(int64_t) * l);
  *qR = malloc(sizeof(int64_t) * l);
  *qinv = malloc(sizeof(int64_t) * l);
}

void rns_gen_q_aux(int64_t* qn, int64_t* qR, int64_t* qinv, const uint64_t* q,
                   uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    qinv[i] = mod_inverse(q[i], (1UL << 63));
    uint64_t R = (1UL << 63) - q[i];
    R += (1UL << 63);
    qR[i] = R % q[i];
    if (qR[i] > R / 2) qR[i] -= q[i];
    qn[i] = mod_inverse(n, q[i]);                                 // 1/n
    qn[i] = mod_mul(qn[i], qR[i], (int64_t)q[i], qinv[i]);        // R/n
    qn[i] = mod_mul_exact(qn[i], qR[i], (int64_t)q[i], qinv[i]);  // R^2/n
  }
}

void rns_gen_q_free(int64_t** qn, int64_t** qR, int64_t** qinv) {
  free((*qinv));
  free((*qR));
  free((*qn));
}
