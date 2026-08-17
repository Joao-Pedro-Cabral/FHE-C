
#include "rns.h"

#include <stdio.h>
#include <stdlib.h>

#include "mod.h"

void rns_new(rns_t* a, const rns_params_t* params) {
  *a = malloc(sizeof(uint64_t*) * params->l);
  for (uint32_t i = 0; i < params->l; i++) {
    poly_new(&(*a)[i], params->n);
  }
}

void rns_free(rns_t* a, const rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_free(&(*a)[i]);
  }
  free(*a);
}

static void rns_gen_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_roots(params->w[i], params->n, (int64_t)params->q[i],
                   params->qR[i], params->qinv[i]);
  }
}

static void rns_gen_inverse_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_inverse_roots(params->w_1[i], params->w[i][params->n / 2],
                           params->n, (int64_t)params->q[i], params->qinv[i]);
  }
}

static inline uint32_t compute_ring_dimension(uint32_t qtotal_bits) {
  if (qtotal_bits < 27) {
    return 1024;
  } else if (qtotal_bits < 54) {
    return 2048;
  } else if (qtotal_bits < 109) {
    return 4096;
  } else if (qtotal_bits < 218) {
    return 8192;
  } else if (qtotal_bits < 438) {
    return 16384;
  } else if (qtotal_bits < 881) {
    return 32768;
  } else if (qtotal_bits < 1747) {
    return 65536;
  } else if (qtotal_bits < 3523) {
    return 131072;
  }
  printf("Too large q: %d\n!", qtotal_bits);
  return 1024;
}

void rns_params_new(rns_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                    uint32_t l) {
  // Allocate Memory
  uint32_t qtotal_bits = q0_bits + l * scale_bits;
  l++;  // q0 prime + l for multiplicative depth
  uint32_t n = compute_ring_dimension(qtotal_bits);
  params->n = n;
  params->l = l;
  params->q = malloc(sizeof(int64_t) * l);
  params->qR2n = malloc(sizeof(int64_t) * l);
  params->qRn = malloc(sizeof(int64_t) * l);
  params->qR = malloc(sizeof(int64_t) * l);
  params->qR2 = malloc(sizeof(int64_t) * l);
  params->qRq = malloc(sizeof(int64_t) * l);
  params->qinv = malloc(sizeof(int64_t) * l);
  params->scales = malloc(sizeof(int64_t) * l);
  rns_new(&params->w, params);
  rns_new(&params->w_1, params);
  // 2n: negacyclic
  find_primes(params->q, q0_bits, scale_bits, 2 * n, l);
  // Init data (less roots)
  for (uint32_t i = 0; i < l; i++) {
    // qinv
    params->qinv[i] = (int64_t)mod64_pow_ct(
        params->q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    // qR
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - params->q[i]);
    params->qR[i] = (int64_t)(R % params->q[i]);
    if (params->qR[i] > (int64_t)(params->q[i] / 2))
      params->qR[i] -= (int64_t)params->q[i];
    // qR2
    params->qR2[i] =
        (int64_t)(((__int128_t)params->qR[i] * params->qR[i]) % params->q[i]);
    if (params->qR2[i] > (int64_t)(params->q[i] / 2))
      params->qR2[i] -= (int64_t)params->q[i];
    if (params->qR2[i] < (int64_t) - (params->q[i] / 2))
      params->qR2[i] += (int64_t)params->q[i];
    // qRn
    params->qRn[i] =
        (int64_t)mod_pow_ct(n, params->q[i] - 2, (int64_t)params->q[i],
                            params->qR[i], params->qinv[i]);  // 1/n
    params->qRn[i] = (int64_t)(((__int128_t)params->qR[i] * params->qRn[i]) %
                               params->q[i]);  // R/n
    if (params->qRn[i] > (int64_t)(params->q[i] / 2))
      params->qRn[i] -= (int64_t)params->q[i];
    if (params->qRn[i] < (int64_t) - (params->q[i] / 2))
      params->qRn[i] += (int64_t)params->q[i];
    // qR2n
    params->qR2n[i] = (int64_t)(((__int128_t)params->qR[i] * params->qRn[i]) %
                                params->q[i]);  // R^2/n
    if (params->qR2n[i] > (int64_t)(params->q[i] / 2))
      params->qR2n[i] -= (int64_t)params->q[i];
    if (params->qR2n[i] < (int64_t) - (params->q[i] / 2))
      params->qR2n[i] += (int64_t)params->q[i];
    // qRq
    params->qRq[i] =
        (int64_t)((((__uint128_t)1 << 64) + params->q[i] / 2) / params->q[i]);
  }
  // scale
  params->scales[l - 1] = (1UL << (scale_bits - 1)) + (1UL << (scale_bits - 2));
  for (uint32_t i = l - 1; i > 0; i--) {
    double scale2 = (double)params->scales[i] * (double)params->scales[i];
    params->scales[i - 1] = (uint64_t)((scale2 + (double)(params->q[i] / 2)) /
                                       (double)params->q[i]);
  }
  // roots
  rns_gen_roots(params);
  rns_gen_inverse_roots(params);
}

void rns_params_free(rns_params_t* params) {
  free(params->qinv);
  free(params->qRq);
  free(params->qR);
  free(params->qR2);
  free(params->q);
  free(params->qR2n);
  free(params->qRn);
  free(params->scales);
  rns_free(&params->w, params);
  rns_free(&params->w_1, params);
}
