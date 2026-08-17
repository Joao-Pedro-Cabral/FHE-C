
#include "ciphertext.h"

#include <stdio.h>

void ciphertext_new(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_new(&ct->ab, params);
  ct->level = params->l;
}

void ciphertext_free(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_free(&ct->ab, params);
}

void ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                          const ciphertext_t* ct1,
                                          const ciphertext_t* ct2,
                                          const switch_t r,
                                          const rns_params_t* params) {
  if (ct1->level != ct2->level) {
    printf(
        "Illegal ciphertext multiplication! Both operands must have same "
        "CRT!\n");
  }
  poly_t tmp;
  rns_t d2;
  poly_new(&tmp, params->n);
  rns_new(&d2, params);
  // Basic Multiplication
  for (uint32_t i = 0; i < ct1->level; i++) {
    // D0
    poly_pointwise_mul(ct3->ab.b[i], ct1->ab.b[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    poly_mulx(ct3->ab.b[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // 1
    // D1
    poly_pointwise_mul(ct3->ab.a[i], ct1->ab.b[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // 1
    poly_pointwise_mul(tmp, ct1->ab.a[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);                        // 1
    poly_add(ct3->ab.a[i], tmp, params->q[i], params->qRq[i], params->n);  // 1
    poly_mulx(ct3->ab.a[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // R
    // D2
    poly_pointwise_mul(d2[i], ct1->ab.a[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // R
  }
  // Relinearization
  for (uint32_t i = 0; i < ct1->level; i++) {
    for (uint32_t j = 0; j < ct1->level; j++) {
      poly_pointwise_mul(tmp, d2[i], r[i].b[j], params->q[j], params->qinv[j],
                         params->n);  // 1
      poly_add(ct3->ab.b[j], tmp, params->q[j], params->qRq[j],
               params->n);  // 1
      poly_pointwise_mul(tmp, d2[i], r[i].a[j], params->q[j], params->qinv[j],
                         params->n);  // R
      poly_add(ct3->ab.a[j], tmp, params->q[j], params->qRq[j],
               params->n);  // R
    }
  }
  // Rescaling (ModSwitch)
  ct3->level = ct1->level - 1;  // We lost one level at each multiplication
  // For ModSwitch we need to come back to the normal domain
  poly_intt(ct3->ab.a[ct3->level], params->w_1[ct3->level],
            params->q[ct3->level], params->qRq[ct3->level],
            params->qinv[ct3->level], params->n);  // R
  poly_mulx(ct3->ab.a[ct3->level], params->qRn[ct3->level],
            params->q[ct3->level], params->qinv[ct3->level],
            params->n);  // R
  for (uint32_t i = 0; i < ct3->level; i++) {
    poly_intt(ct3->ab.a[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);  // R
    poly_mulx(ct3->ab.a[i], params->qRn[i], params->q[i], params->qinv[i],
              params->n);  // R
    poly_mod(tmp, ct3->ab.a[ct3->level], params->q[i], params->qRq[i],
             params->n);
    poly_sub(ct3->ab.a[i], tmp, params->q[i], params->qRq[i], params->n);
    poly_mulx(ct3->ab.a[i], params->qinvq[ct3->level * params->l + i],
              params->q[i], params->qinv[i], params->n);
    poly_ntt(ct3->ab.a[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);  // R
  }
  poly_intt(ct3->ab.b[ct3->level], params->w_1[ct3->level],
            params->q[ct3->level], params->qRq[ct3->level],
            params->qinv[ct3->level], params->n);  // 1
  poly_mulx(ct3->ab.b[ct3->level], params->qRn[ct3->level],
            params->q[ct3->level], params->qinv[ct3->level],
            params->n);  // 1
  for (uint32_t i = 0; i < ct3->level; i++) {
    poly_intt(ct3->ab.b[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);  // R
    poly_mulx(ct3->ab.b[i], params->qRn[i], params->q[i], params->qinv[i],
              params->n);  // 1
    poly_mod(tmp, ct3->ab.b[ct3->level], params->q[i], params->qRq[i],
             params->n);
    poly_sub(ct3->ab.b[i], tmp, params->q[i], params->qRq[i], params->n);
    poly_mulx(ct3->ab.b[i], params->qinvq[ct3->level * params->l + i],
              params->q[i], params->qinv[i], params->n);
    poly_ntt(ct3->ab.b[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);  // 1
  }
  poly_free(&tmp);
  rns_free(&d2, params);
}
