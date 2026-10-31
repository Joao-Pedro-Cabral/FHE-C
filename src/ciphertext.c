
#include "ciphertext.h"

#include <string.h>

#include "mod.h"

void ciphertext_new(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_new(&ct->ab, params);
  ct->level = params->l;
}

void ciphertext_free(ciphertext_t* ct) { rlwe_free(&ct->ab); }

void ciphertext_ciphertext_addition(ciphertext_t* ct3, const ciphertext_t* ct1,
                                    const ciphertext_t* ct2,
                                    const rns_params_t* params) {
  ct3->level = ct1->level;
  uint32_t n = params->n;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_add2(&ct3->ab.a[n * i], &ct2->ab.a[n * i], &ct1->ab.a[n * i],
              params->q[i], params->R_q[i], n);
    poly_add2(&ct3->ab.b[n * i], &ct2->ab.b[n * i], &ct1->ab.b[n * i],
              params->q[i], params->R_q[i], n);
  }
}

void ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                          const ciphertext_t* ct1,
                                          const ciphertext_t* ct2,
                                          const switch_t r,
                                          const rns_params_t* params) {
  poly_t tmp;
  rns_t tmp_r;
  rns_t d2;
  uint32_t n = params->n;
  poly_new(&tmp, n);
  rns_new(&tmp_r, params);
  rns_new(&d2, params);
  // Basic Multiplication
  for (uint32_t i = 0; i < ct1->level; i++) {
    // D0
    poly_pointwise_mul(&ct3->ab.b[n * i], &ct1->ab.b[n * i], &ct2->ab.b[n * i],
                       params->q[i], params->qinvmR[i], n);  // R-1
    // D1
    poly_pointwise_mul(&ct3->ab.a[n * i], &ct1->ab.b[n * i], &ct2->ab.a[n * i],
                       params->q[i], params->qinvmR[i], n);  // 1
    poly_pointwise_mul(tmp, &ct1->ab.a[n * i], &ct2->ab.b[n * i], params->q[i],
                       params->qinvmR[i], n);  // 1
    poly_add(&ct3->ab.a[n * i], tmp, params->q[i], params->R_q[i],
             n);  // 1
    // D2
    poly_pointwise_mul(&d2[n * i], &ct1->ab.a[n * i], &ct2->ab.a[n * i],
                       params->q[i], params->qinvmR[i], n);  // R
    poly_mulx(&d2[n * i], 1, params->q[i], params->qinvmR[i],
              n);  // 1
  }
  // Relinearization
  // We don't need hybrid method for A
  for (uint32_t i = 0; i < ct1->level; i++) {
    for (uint32_t j = 0; j < ct1->level; j++) {
      poly_pointwise_mul(tmp, &d2[n * i], &r[i].a[n * j], params->q[j],
                         params->qinvmR[j],
                         n);  // 1
      poly_add(&ct3->ab.a[n * j], tmp, params->q[j], params->R_q[j],
               n);  // 1
    }
  }
  // Hybrid method for B
  rns_mod_raise(d2, ct1->level, params);
  memset(tmp_r, 0, (params->l + params->k) * params->n * sizeof(int64_t));
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    if (i >= ct1->level && i < params->l) continue;
    for (uint32_t j = 0; j < params->l + params->k; j++) {
      if (j >= ct1->level && j < params->l) continue;
      poly_pointwise_mul(tmp, &d2[n * i], &r[i].b[n * j], params->q[j],
                         params->qinvmR[j],
                         n);  // R-1
      poly_add(&tmp_r[n * j], tmp, params->q[j], params->R_q[j],
               n);  // R-1
    }
  }
  rns_mod_switch(d2, tmp_r, ct1->level, params->l, params->k, params->zRbq,
                 params->yRbq, params);
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_add(&ct3->ab.b[n * i], &tmp_r[n * i], params->q[i], params->R_q[i],
             n);  // R-1
    // Add R factor
    poly_mulx(&ct3->ab.b[n * i], params->R2mq[i], params->q[i],
              params->qinvmR[i],
              n);  // 1
  }
  // Rescaling (ModSwitch)
  ct3->level = ct1->level;
  ciphertext_mod_switch(ct3, tmp_r, params);
  poly_free(&tmp);
  rns_free(&tmp_r);
  rns_free(&d2);
}

void ciphertext_plaintext_addition(ciphertext_t* ct2, const ciphertext_t* ct1,
                                   const poly_t m, const rns_params_t* params) {
  ct2->level = ct1->level;
  uint32_t n = params->n;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_mod(&ct2->ab.b[n * i], m, params->q[i], params->R_q[i], n);
    poly_ntt(&ct2->ab.b[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_add(&ct2->ab.b[n * i], &ct1->ab.b[n * i], params->q[i], params->R_q[i],
             n);
    memcpy(&ct2->ab.a[n * i], &ct1->ab.a[n * i], n * sizeof(int64_t));
  }
}

void ciphertext_plaintext_base_multiplication(ciphertext_t* ct2,
                                              const ciphertext_t* ct1,
                                              const poly_t m,
                                              const rns_params_t* params) {
  ct2->level = ct1->level;
  uint32_t n = params->n;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_mod(&ct2->ab.b[n * i], m, params->q[i], params->R_q[i], n);
    poly_ntt(&ct2->ab.b[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_pointwise_mul(&ct2->ab.a[n * i], &ct1->ab.a[n * i], &ct2->ab.b[n * i],
                       params->q[i], params->qinvmR[i], n);
    poly_pointwise_mul(&ct2->ab.b[n * i], &ct1->ab.b[n * i], &ct2->ab.b[n * i],
                       params->q[i], params->qinvmR[i], n);
    poly_mulx(&ct2->ab.b[n * i], params->R2mq[i], params->q[i],
              params->qinvmR[i], n);
  }
}

// BE CAREFUL: HERE ct->ab.a doesn't have the montgomery factor
void ciphertext_mod_switch(ciphertext_t* ct, rns_t tmp,
                           const rns_params_t* params) {
  ct->level = ct->level - 1;
  // Rescaling (Mod Switch)
  rns_mod_switch(ct->ab.b, tmp, ct->level, ct->level, 1,
                 &params->Rmq[ct->level], params->Rmq, params);
  rns_mod_switch(ct->ab.a, tmp, ct->level, ct->level, 1,
                 &params->Rmq[ct->level], params->Rmq, params);
  // Add R factor
  uint32_t n = params->n;
  for (uint32_t i = 0; i < ct->level; i++) {
    poly_mulx(&ct->ab.a[n * i], params->R2mq[i], params->q[i],
              params->qinvmR[i], n);
  }
}

void ciphertext_plaintext_multiplication(ciphertext_t* ct2,
                                         const ciphertext_t* ct1,
                                         const poly_t m,
                                         const rns_params_t* params) {
  rns_t tmp;
  rns_new(&tmp, params);
  ct2->level = ct1->level;
  ciphertext_plaintext_base_multiplication(ct2, ct1, m, params);
  ciphertext_mod_switch(ct2, tmp, params);
  rns_free(&tmp);
}

// Rotate by 5^h
void ciphertext_rotation(ciphertext_t* ct2, const ciphertext_t* ct1, uint32_t h,
                         rns_t tmp_r, poly_t tmp_p, const switch_t r,
                         const rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < ct1->level; i++) {
    // X -> X^(5^h)
    poly_mulx2(tmp_p, &ct1->ab.a[n * i], params->R_nmq[i], params->q[i],
               params->qinvmR[i], n);
    poly_intt(tmp_p, &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_permute(&tmp_r[n * i], tmp_p, h, n);
    poly_ntt(&tmp_r[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_mulx2(tmp_p, &ct1->ab.b[n * i], params->R_nmq[i], params->q[i],
               params->qinvmR[i], n);
    poly_intt(tmp_p, &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_permute(&ct2->ab.b[n * i], tmp_p, h, n);
    poly_ntt(&ct2->ab.b[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
  // Key Switch
  ct2->level = ct1->level;
  rlwe_key_switch(&ct2->ab, tmp_r, r, tmp_p, ct1->level, params);
}
