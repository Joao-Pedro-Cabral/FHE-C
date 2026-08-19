
#include "ciphertext.h"

#include <stdio.h>
#include <string.h>

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
        "CRT level!\n");
  }
  poly_t tmp;
  rns_t d2;
  poly_new(&tmp, params->n);
  rns_new(&d2, params);
  // Basic Multiplication
  for (uint32_t i = 0; i < ct1->level; i++) {
    // Remove R factor
    poly_mulx(ct1->ab.a[i], 1, params->q[i], params->qinv[i],
              params->n);  // R
    // poly_mulx(ct2->ab.a[i], 1, params->q[i], params->qinv[i],
    //           params->n);  // R
    // D0
    poly_pointwise_mul(ct3->ab.b[i], ct1->ab.b[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    poly_mulx(ct3->ab.b[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // 1
    // D1
    poly_pointwise_mul(ct3->ab.a[i], ct1->ab.b[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    poly_pointwise_mul(tmp, ct1->ab.a[i], ct2->ab.b[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    poly_add(ct3->ab.a[i], tmp, params->q[i], params->qRq[i],
             params->n);  // R-1
    poly_mulx(ct3->ab.a[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // 1
    printf("ct3.a[%d][%d]: %ld, ct3.b[%d][%d]: %ld\n", i, 0, ct3->ab.a[i][0], i,
           0, ct3->ab.b[i][0]);
    // D2
    poly_pointwise_mul(d2[i], ct1->ab.a[i], ct2->ab.a[i], params->q[i],
                       params->qinv[i], params->n);  // R-1
    poly_mulx(d2[i], params->qR2[i], params->q[i], params->qinv[i],
              params->n);  // 1
  }
  // rlwe_t rlwe;
  // rlwe_new(&rlwe, params);
  // for (uint32_t i = 0; i < ct1->level; i++) {
  // for (uint32_t k = 0; k < params->n; k++) {
  //   printf("d2[%d][%d]: %ld\n", i, k, d2[i][k]);
  // }
  // memset(ct3->ab.b[i], 0, params->n * sizeof(int64_t));
  // memset(ct3->ab.a[i], 0, params->n * sizeof(int64_t));
  // memset(rlwe.b[i], 0, params->n*sizeof(int64_t));
  // memset(rlwe.a[i], 0, params->n*sizeof(int64_t));
  // }
  // for (uint32_t i = 0; i < ct1->level; i++) {
  //   for (uint32_t j = 0; j < ct1->level; j++) {
  //     poly_add(rlwe.b[i], r[j].b[i], params->q[i], params->qRq[i],
  //     params->n); poly_add(rlwe.a[i], r[j].a[i], params->q[i],
  //     params->qRq[i], params->n);
  //   }
  // }
  // Relinearization
  for (uint32_t i = 0; i < ct1->level; i++) {
    for (uint32_t j = 0; j < ct1->level; j++) {
      poly_pointwise_mul(tmp, d2[i], r[i].b[j], params->q[j], params->qinv[j],
                         params->n);  // R-1
      poly_mulx(tmp, params->qR2[j], params->q[j], params->qinv[j],
                params->n);  // 1
      poly_add(ct3->ab.b[j], tmp, params->q[j], params->qRq[j],
               params->n);  // 1
      poly_pointwise_mul(tmp, d2[i], r[i].a[j], params->q[j], params->qinv[j],
                         params->n);  // 1
      // poly_mulx(tmp, params->qR2[j], params->q[j], params->qinv[j],
      //           params->n);  // R
      poly_add(ct3->ab.a[j], tmp, params->q[j], params->qRq[j],
               params->n);  // 1
      // printf("ct3.a[%d][%d]: %ld, ct3.b[%d][%d]: %ld\n", j, 0, ct3->ab.a[j][0],
      //        j, 0, ct3->ab.b[j][0]);
    }
  }
  // printf("r[0].a[0][0]: %ld, r[1].a[0][0]: %ld\n", r[0].a[0][0], r[1].a[0][0]);
  // printf("r[0].b[0][0]: %ld, r[1].b[0][0]: %ld\n", r[0].b[0][0], r[1].b[0][0]);
  // for (uint32_t i = 0; i < 1; i++) {
  //   for (uint32_t k = 0; k < 1; k++) {
  //     printf("ct1.a[%d][%d]: %ld, ct1.b[%d][%d]: %ld\n", i, k, ct1->ab.a[i][k],
  //            i, k, ct1->ab.b[i][k]);
  //     printf("ct2.a[%d][%d]: %ld, ct2.b[%d][%d]: %ld\n", i, k, ct2->ab.a[i][k],
  //            i, k, ct2->ab.b[i][k]);
  //     printf("ct3.a[%d][%d]: %ld, ct3.b[%d][%d]: %ld\n", i, k, ct3->ab.a[i][k],
  //            i, k, ct3->ab.b[i][k]);
  //   }
  // }
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
  // rlwe_free(&rlwe, params);
}
