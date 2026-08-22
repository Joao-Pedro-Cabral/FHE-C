
#include "ciphertext.h"

#include <string.h>

#include "mod.h"

void ciphertext_new(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_new(&ct->ab, params);
  ct->level = params->l;
}

void ciphertext_free(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_free(&ct->ab, params);
}

void ciphertext_ciphertext_addition(ciphertext_t* ct3, const ciphertext_t* ct1,
                                    const ciphertext_t* ct2,
                                    const rns_params_t* params) {
  ct3->level = ct1->level;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_add2(ct3->ab.a[i], ct2->ab.a[i], ct1->ab.a[i], params->q[i],
              params->qRq[i], params->n);
    poly_add2(ct3->ab.b[i], ct2->ab.b[i], ct1->ab.b[i], params->q[i],
              params->qRq[i], params->n);
  }
}

void ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                          const ciphertext_t* ct1,
                                          const ciphertext_t* ct2,
                                          const switch_t r,
                                          const rns_params_t* params) {
  poly_t tmp;
  rns_t d2;
  poly_new(&tmp, params->n);
  rns_new(&d2, params);
  // Basic Multiplication
  for (uint32_t i = 0; i < ct1->level; i++) {
    // D0
    poly_pointwise_mul(ct3->ab.b[i], ct1->ab.b[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    // D1
    poly_pointwise_mul(ct3->ab.a[i], ct1->ab.b[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // 1
    poly_pointwise_mul(tmp, ct1->ab.a[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);  // 1
    poly_add(ct3->ab.a[i], tmp, params->q[i], params->qRq[i],
             params->n);  // 1
    // D2
    poly_pointwise_mul(d2[i], ct1->ab.a[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // R
    poly_mulx(d2[i], 1, params->q[i], params->qinv[i],
              params->n);  // 1
  }
  // Relinearization
  for (uint32_t i = 0; i < ct1->level; i++) {
    for (uint32_t j = 0; j < ct1->level; j++) {
      poly_pointwise_mul(tmp, d2[i], r[i].b[j], params->q[j], params->qinv[j],
                         params->n);  // R-1
      poly_add(ct3->ab.b[j], tmp, params->q[j], params->qRq[j],
               params->n);  // R-1
      poly_pointwise_mul(tmp, d2[i], r[i].a[j], params->q[j], params->qinv[j],
                         params->n);  // 1
      poly_add(ct3->ab.a[j], tmp, params->q[j], params->qRq[j],
               params->n);  // 1
    }
  }
  for (uint32_t i = 0; i < ct1->level; i++) {
    // Add R factor
    poly_mulx(ct3->ab.b[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // 1
  }
  // Rescaling (ModSwitch)
  ct3->level = ct1->level - 1;  // We lost one level at each multiplication
  rns_mod_switch(ct3->ab.a, tmp, ct3->level, params);
  rns_mod_switch(ct3->ab.b, tmp, ct3->level, params);
  for (uint32_t i = 0; i < ct1->level; i++) {
    // Add R factor
    poly_mulx(ct3->ab.a[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // R
  }
  poly_free(&tmp);
  rns_free(&d2, params);
}

void ciphertext_plaintext_addition(ciphertext_t* ct2, const ciphertext_t* ct1,
                                   const poly_t m, const rns_params_t* params) {
  ct2->level = ct1->level;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_mod(ct2->ab.b[i], m, params->q[i], params->qRq[i], params->n);
    poly_ntt(ct2->ab.b[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
    poly_add(ct2->ab.b[i], ct1->ab.b[i], params->q[i], params->qRq[i],
             params->n);
    memcpy(ct2->ab.a[i], ct1->ab.a[i], params->n * sizeof(int64_t));
  }
}

void ciphertext_plaintext_multiplication(ciphertext_t* ct2,
                                         const ciphertext_t* ct1,
                                         const poly_t m,
                                         const rns_params_t* params) {
  poly_t tmp;
  poly_new(&tmp, params->n);
  ct2->level = ct1->level - 1;
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_mod(ct2->ab.b[i], m, params->q[i], params->qRq[i], params->n);
    poly_ntt(ct2->ab.b[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
    poly_pointwise_mul(ct2->ab.a[i], ct1->ab.a[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);
    poly_pointwise_mul(ct2->ab.b[i], ct1->ab.b[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);
    poly_mulx(ct2->ab.b[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);
  }
  // Rescaling (Mod Switch)
  rns_mod_switch(ct2->ab.b, tmp, ct2->level, params);
  rns_mod_switch(ct2->ab.a, tmp, ct2->level, params);
  for (uint32_t i = 0; i < ct1->level; i++) {
    // Add R factor
    poly_mulx(ct2->ab.a[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);
  }
  poly_free(&tmp);
}

// Rotate by 5^h
void ciphertext_rotation(ciphertext_t* ct2, ciphertext_t* ct1, uint32_t h,
                         rns_t tmp_r, poly_t tmp_p, const switch_t r,
                         const rns_params_t* params) {
  uint32_t jh = mod_pow2(5, h, 2 * params->n);
  for (uint32_t i = 0; i < ct1->level; i++) {
    // X -> X^(5^h)
    poly_intt(ct1->ab.a[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_mulx(ct1->ab.a[i], params->qRn[i], params->q[i], params->qinv[i],
              params->n);
    poly_permute(tmp_r[i], ct1->ab.a[i], jh, params->n);
    poly_ntt(tmp_r[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
    poly_intt(ct1->ab.b[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_mulx(ct1->ab.b[i], params->qRn[i], params->q[i], params->qinv[i],
              params->n);
    poly_permute(ct2->ab.b[i], ct1->ab.b[i], jh, params->n);
    poly_ntt(ct2->ab.b[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
  }
  // Key Switch
  ct2->level = ct1->level;
  rlwe_key_switch(&ct2->ab, tmp_r, r, tmp_p, ct1->level, params);
}
