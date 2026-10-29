
#include "rns.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mod.h"

void rns_new(rns_t* a, const rns_params_t* params) {
  *a = malloc(sizeof(uint64_t*) * (params->l + params->k) * params->n);
}

void rns_free(rns_t* a) { free(*a); }

// x is in the normal domain
static void rns_fastb_conv(rns_t o, const rns_t x, const uint64_t* q,
                           const uint64_t* b, const int64_t* z,
                           const int64_t* y, const int64_t* qinvmR,
                           const int64_t* binvmR, const int64_t* R_b,
                           uint32_t n, uint32_t k, uint32_t l) {
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      int64_t acc = 0;
      for (uint32_t ii = 0; ii < k; ii++) {
        int64_t tmp = x[n * ii + j];
        tmp = mod_mul_ct(tmp, z[ii], (int64_t)q[ii], qinvmR[ii]);
        tmp = mod_mul_ct(tmp, y[k * i + ii], (int64_t)b[i], binvmR[i]);
        acc = mod_ct(acc + tmp, (int64_t)b[i], R_b[i]);
      }
      o[n * i + j] = acc;
    }
  }
}

void rns_mod_switch(rns_t a, rns_t tmp, uint32_t l, uint32_t m, uint32_t k,
                    const int64_t* z, const int64_t* y,
                    const rns_params_t* params) {
  // For ModSwitch we need to come back to the normal domain
  uint32_t n = params->n;
  // To normal domain
  for (uint32_t i = m; i < m + k; i++) {
    poly_intt(&a[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&a[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
  }
  // FastBConv
  rns_fastb_conv(tmp, &a[n * m], &params->q[m], params->q, z, y,
                 &params->qinvmR[m], params->qinvmR, params->R_q, params->n, k,
                 l);
  for (uint32_t i = 0; i < l; i++) {
    poly_intt(&a[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&a[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    poly_sub(&a[n * i], &tmp[n * i], params->q[i], params->R_q[i], n);
    for (uint32_t j = 0; j < k; j++) {
      poly_mulx(&a[n * i],
                params->Rqinvmq[(m + j) * (params->l + params->k) + i],
                params->q[i], params->qinvmR[i], n);
    }
    poly_ntt(&a[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
}

void rns_mod_raise(rns_t a, uint32_t level, const rns_params_t* params) {
  rns_fastb_conv(&a[params->n * params->l], a, params->q, &params->q[params->l],
                 &params->zRqb[params->l * (level - 1)],
                 &params->yRqb[params->l * (level - 1) * params->k],
                 params->qinvmR, &params->qinvmR[params->l],
                 &params->R_q[params->l], params->n, level, params->k);
}

static void rns_gen_roots(rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    poly_gen_roots(&params->w[n * i], n, (int64_t)params->q[i], params->Rmq[i],
                   params->qinvmR[i]);
  }
}

static void rns_gen_inverse_roots(rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    poly_gen_inverse_roots(&params->w_1[n * i], params->w[n * i + n / 2], n,
                           (int64_t)params->q[i], params->qinvmR[i]);
  }
}

static inline uint32_t compute_ring_dimension(uint32_t qtotal_bits) {
  return 128;
}

void rns_params_new(rns_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                    uint32_t l) {
  // Allocate Memory
  uint32_t qtotal_bits = q0_bits + l * scale_bits;
  l++;  // q0 prime + l for multiplicative depth
  uint32_t k = P_MOD_SIZE;
  uint32_t n = compute_ring_dimension(qtotal_bits);
  params->n = n;
  params->l = l;
  params->k = P_MOD_SIZE;
  uint32_t lk = l + P_MOD_SIZE;
  params->q = malloc(sizeof(int64_t) * lk);
  params->R2_nmq = malloc(sizeof(int64_t) * lk);
  params->_nmq = malloc(sizeof(int64_t) * lk);
  params->R_nmq = malloc(sizeof(int64_t) * lk);
  params->Rmq = malloc(sizeof(int64_t) * lk * k);  // k times for ModSwitch
  params->R2mq = malloc(sizeof(int64_t) * lk);
  params->R_q = malloc(sizeof(int64_t) * lk);
  params->qinvmR = malloc(sizeof(int64_t) * lk);
  params->Rqinvmq = malloc(sizeof(int64_t) * lk * lk);  // only for rescaling
  params->scales = malloc(sizeof(int64_t) * l);  // Scales only for Q moduli
  params->zRbq = malloc(sizeof(int64_t) * k);
  params->yRbq = malloc(sizeof(int64_t) * k * l);
  params->zRqb = malloc(sizeof(int64_t) * l * l);
  params->yRqb = malloc(sizeof(int64_t) * k * l * l);
  rns_new(&params->w, params);
  rns_new(&params->w_1, params);
  // 2n: negacyclic
  find_primes(params->q, q0_bits, scale_bits, P_MOD_BITS, 2 * n, l, P_MOD_SIZE);
  // Init data (less roots)
  for (uint32_t i = 0; i < lk; i++) {
    // qinvmR
    params->qinvmR[i] = (int64_t)mod64_pow_ct(
        params->q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    // Rmq
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - params->q[i]);
    params->Rmq[i] =
        q_to_q2((int64_t)(R % params->q[i]), (int64_t)(params->q[i]));
    // R2mq
    params->R2mq[i] = q_to_q2(
        (int64_t)(((__int128_t)params->Rmq[i] * params->Rmq[i]) % params->q[i]),
        (int64_t)(params->q[i]));
    // _nmq
    params->_nmq[i] =
        q_to_q2(mod_pow_ct(n, params->q[i] - 2, (int64_t)params->q[i],
                           params->Rmq[i], params->qinvmR[i]),
                (int64_t)(params->q[i]));  // 1/n
    // R_nmq
    params->R_nmq[i] =
        q_to_q2((int64_t)(((__int128_t)params->Rmq[i] * params->_nmq[i]) %
                          params->q[i]),
                (int64_t)(params->q[i]));  // R/n
    // R2_nmq
    params->R2_nmq[i] =
        q_to_q2((int64_t)(((__int128_t)params->Rmq[i] * params->R_nmq[i]) %
                          params->q[i]),
                (int64_t)(params->q[i]));  // R^2/n
    // R_q
    params->R_q[i] =
        (int64_t)((((__uint128_t)1 << 64) + params->q[i] / 2) / params->q[i]);
  }
  // Complete Rmq
  for (uint32_t i = 1; i < k; i++) {
    memcpy(&params->Rmq[i * lk], params->Rmq, sizeof(int64_t) * lk);
  }
  // Rqinvmq
  for (uint32_t j = 0; j < lk; j++) {
    for (uint32_t i = 0; i < lk; i++) {
      if (i == j) {
        params->Rqinvmq[i * lk + j] = 0;  // placeholder (never used)
      } else {
        params->Rqinvmq[i * lk + j] =
            mod_pow_ct((int64_t)params->q[i], params->q[j] - 2,
                       (int64_t)params->q[j], params->Rmq[j],
                       params->qinvmR[j]);  // (q_i)-1
        params->Rqinvmq[i * lk + j] =
            q_to_q2((int64_t)(((__int128_t)params->Rmq[j] *
                               params->Rqinvmq[i * lk + j]) %
                              params->q[j]),
                    (int64_t)(params->q[j]));  // Add R factor
      }
    }
  }
  // scale
  params->scales[l - 1] = (1UL << (scale_bits - 1)) + (1UL << (scale_bits - 2));
  for (uint32_t i = l - 1; i > 0; i--) {
    double scale2 = (double)params->scales[i] * (double)params->scales[i];
    params->scales[i - 1] = (uint64_t)((scale2 + (double)(params->q[i] / 2)) /
                                       (double)params->q[i]);
  }
  // zRbq
  for (uint32_t i = params->l; i < lk; i++) {
    int64_t acc = 1;
    for (uint32_t j = params->l; j < lk; j++) {
      if (i == j) continue;
      acc =
          mod_mul_ct(mod_mul_ct(acc, (int64_t)params->q[j],
                                (int64_t)params->q[i], params->qinvmR[i]),
                     params->R2mq[i], (int64_t)params->q[i], params->qinvmR[i]);
    }
    acc = mod_pow_ct(acc, params->q[i] - 2, (int64_t)params->q[i],
                     params->Rmq[i], params->qinvmR[i]);
    acc = mod_mul_ct(acc, params->R2mq[i], (int64_t)params->q[i],
                     params->qinvmR[i]);
    params->zRbq[i - params->l] = q_to_q2(acc, (int64_t)(params->q[i]));
  }
  // yRbq
  for (uint32_t i = 0; i < params->l; i++) {
    for (uint32_t j = params->l; j < lk; j++) {
      int64_t acc = 1;
      for (uint32_t ii = params->l; ii < lk; ii++) {
        if (j == ii) continue;
        acc = mod_mul_ct(mod_mul_ct(acc, (int64_t)params->q[ii],
                                    (int64_t)params->q[i], params->qinvmR[i]),
                         params->R2mq[i], (int64_t)params->q[i],
                         params->qinvmR[i]);
      }
      params->yRbq[params->k * i + j - params->l] =
          q_to_q2(acc, (int64_t)(params->q[i]));
    }
  }
  // zRqb
  // Base cases
  params->zRqb[0] = params->Rmq[0];  // Not used
  params->zRqb[params->l + 0] = params->Rqinvmq[lk];
  params->zRqb[params->l + 1] = params->Rqinvmq[1];
  for (uint32_t i = 2; i < params->l; i++) {
    for (uint32_t j = 0; j < i; j++) {
      int64_t tmp = mod_mul_ct(params->zRqb[(i - 1) * params->l + j],
                               params->Rqinvmq[i * lk + j],
                               (int64_t)params->q[j], params->qinvmR[j]);
      params->zRqb[i * params->l + j] = q_to_q2(tmp, (int64_t)(params->q[i]));
    }
    int64_t acc = params->Rqinvmq[i];
    for (uint32_t j = 1; j < i; j++) {
      acc = mod_mul_ct(acc, params->Rqinvmq[j * lk + i], (int64_t)params->q[i],
                       params->qinvmR[i]);
    }
    params->zRqb[i * params->l + i] = q_to_q2(acc, (int64_t)(params->q[i]));
  }
  // yRqb
  // Base case
  for (uint32_t i = 0; i < params->k; i++) {
    params->yRqb[i * params->l] = params->Rmq[params->l + i];  // Not used
  }
  for (uint32_t i = 0; i < params->k; i++) {
    uint32_t li = params->l + i;
    for (uint32_t j = 1; j < params->l; j++) {
      int64_t tmp;
      for (uint32_t ii = 0; ii < j; ii++) {
        tmp = mod_mul_ct(
            params->yRqb[(j - 1) * params->l * params->k + i * params->l + ii],
            (int64_t)params->q[j], (int64_t)params->q[li], params->qinvmR[li]);
        tmp = mod_mul_ct(tmp, params->R2mq[li], (int64_t)params->q[li],
                         params->qinvmR[li]);
        params->yRqb[j * params->l * params->k + i * params->l + ii] =
            q_to_q2(tmp, (int64_t)(params->q[li]));
      }
      tmp = mod_mul_ct(
          params
              ->yRqb[(j - 1) * params->l * params->k + i * params->l + (j - 1)],
          (int64_t)params->q[j - 1], (int64_t)params->q[li],
          params->qinvmR[li]);
      tmp = mod_mul_ct(tmp, params->R2mq[li], (int64_t)params->q[li],
                       params->qinvmR[li]);
      params->yRqb[j * params->l * params->k + i * params->l + j] =
          q_to_q2(tmp, (int64_t)(params->q[li]));
    }
  }
  // roots
  rns_gen_roots(params);
  rns_gen_inverse_roots(params);
}

void rns_params_free(rns_params_t* params) {
  free(params->qinvmR);
  free(params->Rqinvmq);
  free(params->R_q);
  free(params->Rmq);
  free(params->R2mq);
  free(params->q);
  free(params->R2_nmq);
  free(params->R_nmq);
  free(params->_nmq);
  free(params->zRbq);
  free(params->yRbq);
  free(params->zRqb);
  free(params->yRqb);
  free(params->scales);
  rns_free(&params->w);
  rns_free(&params->w_1);
}
