
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

void rns_gen_roots(rns roots, const uint64_t* q, const int64_t* qR,
                   const int64_t* qinv, uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_roots(roots[i], n, (int64_t)q[i], qR[i], qinv[i]);
  }
}

void rns_gen_inverse_roots(rns roots, const uint64_t* q, const int64_t* qinv,
                           uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    poly_gen_inverse_roots(roots[i], n, (int64_t)q[i], qinv[i]);
  }
}

void rns_gen_q_new(int64_t** qR2n, int64_t** qR, int64_t** qRq, int64_t** qinv, uint32_t l) {
  *qR2n = malloc(sizeof(int64_t) * l);
  *qR = malloc(sizeof(int64_t) * l);
  *qRq = malloc(sizeof(int64_t) * l);
  *qinv = malloc(sizeof(int64_t) * l);
}

void rns_gen_q_aux(int64_t* qR2n, int64_t* qR, int64_t* qRq, int64_t* qinv, const uint64_t* q,
                   uint32_t n, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    qinv[i] = (int64_t)mod_pow(1UL << 32, 2*q[i]-4, q[i]); // (1/2^64) = (2^64)^(q-2) = (2^32)^(2q-4)
    if (qinv[i] > (int64_t)(q[i] / 2)) qinv[i] -= (int64_t)q[i];
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - q[i]);
    qR[i] = (int64_t)(R % q[i]);
    if (qR[i] > (int64_t)(q[i] / 2)) qR[i] -= (int64_t)q[i];
    qR2n[i] = (int64_t)mod_pow(n, q[i]-2, q[i]);                             // 1/n
    qR2n[i] = (int64_t)((((__int128_t)1 << 64) * qR2n[i]) % q[i]);  // R/n
    qR2n[i] = (int64_t)((((__int128_t)1 << 64) * qR2n[i]) % q[i]);  // R^2/n
    if (qR2n[i] > (int64_t)(q[i] / 2)) qR2n[i] -= (int64_t)q[i];
    if (qR2n[i] < (int64_t) - (q[i] / 2)) qR2n[i] += (int64_t)q[i];
    qRq[i] = (int64_t)((((__uint128_t)1 << 64) + q[i]/2)/q[i]);
  }
}

void rns_gen_q_free(int64_t** qR2n, int64_t** qR, int64_t** qRq, int64_t** qinv) {
  free((*qinv));
  free((*qRq));
  free((*qR));
  free((*qR2n));
}
