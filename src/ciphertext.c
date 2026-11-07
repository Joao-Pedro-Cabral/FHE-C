
#include "ciphertext.h"

#include <assert.h>
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
  const uint32_t n = params->n;
  const uint32_t level = ct1->level;
  const uint32_t total = params->l + params->k;

  assert(level == ct2->level);
  assert(level > 1 && level <= params->l);
  assert(ct3 != ct1 && ct3 != ct2);

  poly_t tmp, tmp2;
  rns_t d2, acc_a, acc_b, scratch;

  poly_new(&tmp, n);
  poly_new(&tmp2, n);

  rns_new(&d2, params);
  rns_new(&acc_a, params);
  rns_new(&acc_b, params);
  rns_new(&scratch, params);

  // ----------------------------------------
  // 1. Basic ciphertext multiplication
  // ----------------------------------------
  // Inputs: b ordinary NTT, a Montgomery NTT
  // Outputs: D0, D1, D2 ordinary NTT

  for (uint32_t i = 0; i < level; i++) {
    uint64_t q = params->q[i];
    int64_t qinv = params->qinvmR[i];

    // D0 = b1 * b2
    poly_pointwise_mul(&ct3->ab.b[n * i], &ct1->ab.b[n * i], &ct2->ab.b[n * i],
                       q, qinv, n);
    poly_mulx(&ct3->ab.b[n * i], params->R2mq[i], q, qinv, n);

    // D1 = b1 * a2 + a1 * b2
    poly_pointwise_mul(&ct3->ab.a[n * i], &ct1->ab.b[n * i], &ct2->ab.a[n * i],
                       q, qinv, n);
    poly_pointwise_mul(tmp, &ct1->ab.a[n * i], &ct2->ab.b[n * i], q, qinv, n);

    poly_add(&ct3->ab.a[n * i], tmp, q, params->R_q[i], n);

    // D2 = a1 * a2
    poly_pointwise_mul(&d2[n * i], &ct1->ab.a[n * i], &ct2->ab.a[n * i], q,
                       qinv, n);
    poly_mulx(&d2[n * i], 1, q, qinv, n);
  }

  // ----------------------------------------
  // 2. Prepare RNS decomposition digits
  // ----------------------------------------
  // Convert D2 to coefficient representation
  // and extend the basis from Q to QP.

  rns_mod_raise(d2, level, params);

  // ----------------------------------------
  // 3. Relinearization over QP
  // ----------------------------------------
  // Compute:
  //   (acc_a, acc_b) = sum_i D2_i * r[i]
  //
  // Only active Q primes are decomposition
  // indices; P primes are destination moduli.

  memset(acc_a, 0, total * n * sizeof(int64_t));
  memset(acc_b, 0, total * n * sizeof(int64_t));

  for (uint32_t i = 0; i < total; i++) {
    for (uint32_t j = 0; j < total; j++) {
      // Skip inactive ciphertext primes.
      if (j >= level && j < params->l) continue;

      uint64_t q = params->q[j];
      int64_t qinv = params->qinvmR[j];

      // Lift digit D2_i from q_i to q_j.
      memcpy(tmp2, &d2[n * i], n * sizeof(int64_t));

      // for (uint32_t t = 0; t < n; t++) {
      //   int64_t x = tmp2[t] % (int64_t)q;
      //   if (x < 0) x += (int64_t)q;
      //   tmp2[t] = x;
      // }

      // NTT in the destination modulus.
      poly_ntt(tmp2, &params->w[n * j], q, params->R_q[j], qinv, n);

      // A: ordinary * Montgomery -> ordinary
      poly_pointwise_mul(tmp, tmp2, &r[i].a[n * j], q, qinv, n);

      poly_add(&acc_a[n * j], tmp, q, params->R_q[j], n);

      // B: ordinary * ordinary -> R^-1
      poly_pointwise_mul(tmp, tmp2, &r[i].b[n * j], q, qinv, n);
      poly_add(&acc_b[n * j], tmp, q, params->R_q[j], n);
    }
  }

  // ----------------------------------------
  // 4. Normalize Montgomery representation
  // ----------------------------------------
  // acc_a is ordinary NTT.
  // acc_b has factor R^-1.

  for (uint32_t j = 0; j < total; j++) {
    if (j >= level && j < params->l) continue;

    poly_mulx(&acc_b[n * j], params->R2mq[j], params->q[j], params->qinvmR[j],
              n);
  }

  // ----------------------------------------
  // 5. ModDown from QP to Q
  // ----------------------------------------
  // Both components are reduced by P.

  rns_mod_switch(acc_a, scratch, level, params->l, params->k, params->zRbq,
                 params->yRbq, params);

  rns_mod_switch(acc_b, scratch, level, params->l, params->k, params->zRbq,
                 params->yRbq, params);

  // ----------------------------------------
  // 6. Complete relinearization
  // ----------------------------------------
  // C0 = D0 + acc_b
  // C1 = D1 + acc_a

  for (uint32_t j = 0; j < level; j++) {
    uint64_t q = params->q[j];

    poly_add(&ct3->ab.b[n * j], &acc_b[n * j], q, params->R_q[j], n);

    poly_add(&ct3->ab.a[n * j], &acc_a[n * j], q, params->R_q[j], n);
  }

  // ----------------------------------------
  // 7. CKKS rescaling
  // ----------------------------------------
  // Drop the last active ciphertext prime.
  // Restore Montgomery factor R to C1.

  ct3->level = level;
  ciphertext_mod_switch(ct3, scratch, params);

  // ----------------------------------------
  // 8. Cleanup
  // ----------------------------------------

  poly_free(&tmp);
  poly_free(&tmp2);

  rns_free(&d2);
  rns_free(&acc_a);
  rns_free(&acc_b);
  rns_free(&scratch);
}

void old_ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                              const ciphertext_t* ct1,
                                              const ciphertext_t* ct2,
                                              const switch_t r,
                                              const rns_params_t* params) {
  poly_t tmp, tmp2;
  rns_t tmp_r;
  rns_t d2;
  uint32_t n = params->n;
  poly_new(&tmp, n);
  poly_new(&tmp2, n);
  rns_new(&tmp_r, params);
  rns_new(&d2, params);
  // Basic Multiplication
  for (uint32_t i = 0; i < ct1->level; i++) {
    // D0
    poly_pointwise_mul(&ct3->ab.b[n * i], &ct1->ab.b[n * i], &ct2->ab.b[n * i],
                       params->q[i], params->qinvmR[i], n);  // R-1
    poly_mulx(&ct3->ab.b[n * i], params->R2mq[i], params->q[i],
              params->qinvmR[i], n);  // 1
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
  // Hybrid method for B
  rns_mod_raise(d2, ct1->level, params);
  memset(tmp_r, 0, (params->l + params->k) * params->n * sizeof(int64_t));
  // We don't need hybrid method for A
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    if (i >= ct1->level && i < params->l) continue;
    for (uint32_t j = 0; j < params->l + params->k; j++) {
      if (j >= ct1->level && j < params->l) continue;
      memcpy(tmp2, &d2[n * i], n * sizeof(int64_t));
      poly_ntt(tmp2, &params->w[n * j], params->q[j], params->R_q[j],
               params->qinvmR[j], n);
      if (i < ct1->level && j < ct1->level) {
        poly_pointwise_mul(tmp, tmp2, &r[i].a[n * j], params->q[j],
                           params->qinvmR[j],
                           n);  // 1
        poly_add(&ct3->ab.a[n * j], tmp, params->q[j], params->R_q[j],
                 n);  // 1
      }
      poly_pointwise_mul(tmp, tmp2, &r[i].b[n * j], params->q[j],
                         params->qinvmR[j], n);  // R-1
      poly_add(&tmp_r[n * j], tmp, params->q[j], params->R_q[j],
               n);  // R-1
    }
  }
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    if (i >= ct1->level && i < params->l) continue;
    // Add R factor
    poly_mulx(&tmp_r[n * i], params->R2mq[i], params->q[i], params->qinvmR[i],
              n);  // 1
  }
  rns_mod_switch(tmp_r, d2, ct1->level, params->l, params->k, params->zRbq,
                 params->yRbq, params);
  for (uint32_t i = 0; i < ct1->level; i++) {
    poly_add(&ct3->ab.b[n * i], &tmp_r[n * i], params->q[i], params->R_q[i],
             n);  // 1
  }
  // Rescaling (ModSwitch)
  ct3->level = ct1->level;
  ciphertext_mod_switch(ct3, tmp_r, params);
  poly_free(&tmp2);
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
                         rns_t tmp_r, rns_t tmp_r2, poly_t tmp_p, poly_t tmp_p2,
                         const switch_t r, const rns_params_t* params) {
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
  rlwe_key_switch(&ct2->ab, tmp_r, r, tmp_p, tmp_p2, tmp_r2, ct1->level,
                  params);
}
