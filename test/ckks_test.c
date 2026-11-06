
#include "ckks.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "complex.h"
#include "mod.h"
#include "random.h"
#include "utils.h"

#ifdef __RUN_MODE__
#define N_TESTS 1
#elif __VTUNE_MODE__
#define N_TESTS 1
#else
#define N_TESTS 1
#endif

#ifdef __RUN_MODE__
static uint64_t time_ns(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}
#endif

void check_q_aux(const uint64_t* q, const int64_t* Rmq, const int64_t* R2mq,
                 const int64_t* R_q, const int64_t* qinvmR, const int64_t* _nmq,
                 const int64_t* R_nmq, const int64_t* R2_nmq,
                 const int64_t* Rqinvmq, const int64_t* zRbq,
                 const int64_t* yRbq, const int64_t* zRqb, const int64_t* yRqb,
                 uint32_t l, uint32_t k) {
  printf("---------- VERIFY Q AUXILIARIES ROOTS ----------\n");
  for (uint32_t i = 0; i < l + k; i++) {
    printf(
        "Rmq: %ld, R2mq: %ld, qinvmR: %ld, R_q: %ld, _nmq: %ld, R_nmq: %ld , "
        "R2_nmq: %ld "
        "mod "
        "%ld\n",
        Rmq[i], R2mq[i], qinvmR[i], R_q[i], _nmq[i], R_nmq[i], R2_nmq[i], q[i]);
  }
  for (uint32_t i = 0; i < l + k; i++) {
    for (uint32_t j = 0; j < l + k; j++) {
      printf("Rqinvmq[%d]: %ld mod %ld\n", (l + k) * i + j,
             Rqinvmq[(l + k) * i + j], q[j]);
    }
  }
  for (uint32_t i = 0; i < k; i++) {
    printf("zRbq[%d]: %ld mod %ld\n", i, zRbq[i], q[l + i]);
  }
  for (uint32_t i = 0; i < k; i++) {
    for (uint32_t j = 0; j < l; j++) {
      printf("yRbq[%d]: %ld mod %ld\n", l * i + j, yRbq[l * i + j], q[j]);
    }
  }
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j <= i; j++) {
      printf("zRqb[%d]: %ld mod %ld\n", l * i + j, zRqb[l * i + j], q[j]);
    }
  }
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < k; j++) {
      for (uint32_t m = 0; m <= i; m++) {
        printf("yRqb[%d]: %ld mod %ld\n", l * k * i + l * j + m,
               yRqb[l * k * i + l * j + m], q[l + j]);
      }
    }
  }
}

void check_roots(const rns_t roots, const uint64_t* q, uint32_t lk,
                 uint32_t n) {
  printf("---------- VERIFY RNS ROOTS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < lk; i++) {
    for (uint32_t j = 0; j < n; j++) {
      if ((roots[n * i + j] >= (int64_t)q[i]) |
          (roots[n * i + j] <= -(int64_t)q[i])) {
        correct = 0;
        printf("roots[%d]: %ld mod %ld\n", j, roots[n * i + j], q[i]);
      }
    }
  }
  if (correct) {
    printf("---------- RNS ROOTS CORRECT ----------\n");
  } else {
    printf("---------- RNS ROOTS WRONG ----------\n");
  }
}

void check_secret_key(const ckks_params_t* params) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly_t s_hat;
  uint32_t n = params->rns.n;
  poly_new(&s_hat, n);
  for (uint32_t i = 0; i < params->rns.l; i++) {
    memcpy(s_hat, &params->keys.s[n * i], n * sizeof(uint64_t));
    poly_intt(s_hat, &params->rns.w_1[n * i], params->rns.q[i],
              params->rns.R_q[i], params->rns.qinvmR[i], n);
    poly_mulx(s_hat, params->rns.R_nmq[i], params->rns.q[i],
              params->rns.qinvmR[i],
              n);  // multiply by R^2/n
    for (uint32_t j = 0; j < n; j++) {
      if (s_hat[j] != params->keys.s_hat[j]) {
        correct = 0;
        printf("Invalid secret key! s[%d]: %ld instead of %ld mod %ld\n", j,
               s_hat[j], params->keys.s_hat[j], params->rns.q[i]);
      }
    }
  }
  poly_free(&s_hat);
  if (correct) {
    printf("---------- SECRET KEY CORRECT ----------\n");
  } else {
    printf("---------- SECRET KEY WRONG ----------\n");
  }
}

void check_relinearization_key(const ckks_params_t* params) {
  printf("---------- VERIFY RELINEARIZATION KEY ----------\n");
  uint8_t correct = 1;
  poly_t s2_hat, s2, s2_aux, e;
  poly_new(&s2_hat, 2 * params->rns.n);
  poly_new(&s2, params->rns.n);
  poly_new(&s2_aux, params->rns.n);
  poly_new(&e, params->rns.n);
  memset(s2_hat, 0, 2 * params->rns.n * sizeof(int64_t));
  // S^2
  // Since the coefficients are ternary and q_0 >> N ( > 2*N)
  // Then the products/sums will never wrap around on any moduli
  for (size_t i = 0; i < params->rns.n; i++) {
    for (size_t j = 0; j < params->rns.n; j++) {
      s2_hat[i + j] += params->keys.s_hat[i] * params->keys.s_hat[j];
    }
  }
  // Reduction over X^n + 1
  for (size_t i = 0; i < params->rns.n; i++) {
    s2_hat[i] = s2_hat[i] - s2_hat[params->rns.n + i];
  }
  rlwe_key_switch_decrypt(s2, params->keys.r, &params->keys, &params->rns);
  for (uint32_t j = 0; j < params->rns.n; j++) {
    if (s2[j] != s2_hat[j]) {
      correct = 0;
      printf("Invalid reli. key! s2[%d]: %ld, s2_hat[%d]: %ld mod %ld\n", j,
             s2[j], j, s2_hat[j], params->rns.q[0]);
    }
  }
  if (correct) {
    printf("---------- RELINEARIZATION KEY CORRECT ----------\n");
  } else {
    printf("---------- RELINEARIZATION KEY WRONG ----------\n");
  }
  poly_free(&s2_hat);
  poly_free(&s2_aux);
  poly_free(&s2);
  poly_free(&e);
}

void check_rotation_key(const ckks_params_t* params) {
  printf("---------- VERIFY ROTATION KEY ----------\n");
  uint8_t correct = 1;
  poly_t s2;
  poly_new(&s2, params->rns.n);
  // Since the rotation multiplies by +-1 the number of -+1 and 0 is constant
  uint32_t cnt0 = 0, cnt1 = 0;
  uint32_t cnt20 = 0, cnt21 = 0;
  for (size_t i = 0; i < params->rns.n; i++) {
    if (params->keys.s_hat[i] == 0)
      cnt0++;
    else
      cnt1++;
  }
  uint32_t log2_n = floor_log2(params->rns.n) - 1;
  for (uint32_t i = 0; i < log2_n; i++) {
    cnt20 = cnt21 = 0;
    rlwe_key_switch_decrypt(s2, params->keys.ip[i], &params->keys,
                            &params->rns);
    for (uint32_t k = 0; k < params->rns.n; k++) {
      if (s2[k] == 0)
        cnt20++;
      else
        cnt21++;
    }
    if ((cnt20 != cnt0) || (cnt21 != cnt1)) {
      correct = 0;
      printf(
          "Invalid rot. key count! cnt0: %d, cnt20: %d, cnt1: %d, cnt21: %d, "
          "for %d left rotation mod %ld\n",
          cnt0, cnt20, cnt1, cnt21, 1 << i, params->rns.q[i]);
    }
  }
  if (correct) {
    printf("---------- ROTATION KEY CORRECT ----------\n");
  } else {
    printf("---------- ROTATION KEY WRONG ----------\n");
  }
  poly_free(&s2);
}

void check_keys(const ckks_params_t* params) {
  check_secret_key(params);
  // check_relinearization_key(params);
  check_rotation_key(params);
}

void check_ciphertext(const ciphertext_t* ct, const rns_params_t* params) {
  printf("---------- VERIFY CYPHERTEXT ----------\n");
  uint8_t correct = 1;
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      if ((ct->ab.a[n * i + j] > (int64_t)params->q[i] / 2) |
          (ct->ab.a[n * i + j] < -(int64_t)params->q[i] / 2)) {
        correct = 0;
        printf("a[%d][%d]: %ld mod %ld\n", i, j, ct->ab.a[n * i + j],
               params->q[i]);
      }
      if ((ct->ab.b[n * i + j] > (int64_t)params->q[i] / 2) |
          (ct->ab.b[n * i + j] < -(int64_t)params->q[i] / 2)) {
        correct = 0;
        printf("b[%d][%d]: %ld mod %ld\n", i, j, ct->ab.b[n * i + j],
               params->q[i]);
      }
    }
  }
  if (correct) {
    printf("---------- CYPHERTEXT CORRECT ----------\n");
  } else {
    printf("---------- CYPHERTEXT WRONG ----------\n");
  }
}

void check_encoded_polynomial(const poly_t m, const int64_t* q, uint32_t n) {
  printf("---------- VERIFY ENCODED POLYNOMIAL ----------\n");
  uint8_t correct = 1;
  int64_t q0 = q[0];
  for (uint32_t i = 0; i < n; i++) {
    if ((m[i] > (int64_t)q0 / 2) | (m[i] < -(int64_t)q0 / 2)) {
      correct = 0;
      printf("m[%d]: %ld mod %ld\n", i, m[i], q0);
    }
  }
  if (correct) {
    printf("---------- ENCODED POLYNOMIAL CORRECT ----------\n");
  } else {
    printf("---------- ENCODED POLYNOMIAL WRONG ----------\n");
  }
}

void check_decoded_polynomial(const poly_t m, const poly_t m2, uint64_t scale,
                              uint32_t n) {
  printf("---------- VERIFY DECODED POLYNOMIAL ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < n; i++) {
    int64_t diff = m[i] - m2[i];
    if ((diff > (int64_t)scale / 2) | (diff < -(int64_t)scale / 2)) {
      correct = 0;
      printf("m[%d]: %ld, m2[%d]: %ld\n", i, m[i], i, m2[i]);
    }
  }
  if (correct) {
    printf("---------- DECODED POLYNOMIAL CORRECT ----------\n");
  } else {
    printf("---------- DECODED POLYNOMIAL WRONG ----------\n");
  }
}

void check_vectors(const vec_t* v, const vec_t* v2, uint32_t n) {
  printf("---------- VERIFY VECTORS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < n; i++) {
    double diff = v->rel[i] - v2->rel[i];
    // TODO: Better metric
    if ((diff > 0.5) | (diff < -0.5)) {
      correct = 0;
      printf("v_rel[%d]: %lf, v2_rel[%d]: %lf\n", i, v->rel[i], i, v2->rel[i]);
    }
    diff = v->img[i] - v2->img[i];
    // TODO: Better metric
    if ((diff > 0.5) | (diff < -0.5)) {
      correct = 0;
      printf("v_img[%d]: %lf, v2_img[%d]: %lf\n", i, v->img[i], i, v2->img[i]);
    }
  }
  if (correct) {
    printf("---------- VECTORS CORRECT ----------\n");
  } else {
    printf("---------- VECTORS WRONG ----------\n");
  }
}

bool check_ntt_mod_switch(rns_t a, rns_t a2, const rns_params_t* params) {
  bool correct = true;
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l - 1; i++) {
    poly_intt(&a[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&a[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    for (uint32_t j = 0; j < n; j++) {
      int64_t new_a =
          (int64_t)(((__int128_t)params->q[params->l - 1] * a[n * i + j] +
                     a[(params->l - 1) * n + j]) %
                    (int64_t)params->q[i]);
      if (new_a < 0) new_a += (int64_t)params->q[i];
      int64_t diff = new_a - a2[n * i + j];
      int64_t error_threshold = 2;  // Only removes one moduli
      if ((diff >= error_threshold) || (diff <= -error_threshold)) {
        correct = false;
        printf("Incorrect new_a: %ld, a[%d][%d]: %ld\n", new_a, i, j,
               a2[n * i + j]);
      }
    }
  }
  return correct;
}

void check_mod_switch(const rns_params_t* params) {
  printf("---------- VERIFY MOD SWITCH ----------\n");
  uint8_t correct = 1;
  rns_t a, a2;
  poly_t tmp;
  uint32_t n = params->n;
  rns_new(&a, params);
  rns_new(&a2, params);
  poly_new(&tmp, n);
  for (uint32_t i = 0; i < params->l; i++) {
    randombytes((uint8_t*)a[i], n * sizeof(int64_t));
    for (uint32_t j = 0; j < n; j++) {
      a[n * i + j] = a[n * i + j] % (int64_t)params->q[i];
      if (a[n * i + j] < 0) a[n * i + j] += (int64_t)params->q[i];
    }
  }
  memcpy(a2, a, n * params->l * sizeof(int64_t));
  for (uint32_t i = 0; i < params->l; i++) {
    poly_ntt(&a[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
  rns_mod_switch(a, tmp, params->l - 1, params->l - 1, 1,
                 &params->Rmq[params->l - 1], params->Rmq, params);
  correct = check_ntt_mod_switch(a, a2, params);
  poly_free(&tmp);
  rns_free(&a);
  rns_free(&a2);
  if (correct) {
    printf("---------- MOD SWITCH CORRECT ----------\n");
  } else {
    printf("---------- MOD SWITCH WRONG ----------\n");
  }
}

void check_d0(const ciphertext_t* ct2, const ciphertext_t* ct,
              const rns_params_t* params) {
  printf("---------- VERIFY D0 ----------\n");
  uint8_t correct = 1;

  rns_t d0;
  rns_new(&d0, params);
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l; i++) {
    // D0 = B1^2
    for (uint32_t j = 0; j < n; j++) {
      d0[n * i + j] =
          (int64_t)(((__int128_t)ct->ab.b[n * i + j] * ct->ab.b[n * i + j]) %
                    (int64_t)params->q[i]);
    }
    poly_intt(&d0[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&d0[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    for (uint32_t j = 0; j < n; j++) {
      if (d0[n * i + j] < 0) d0[n * i + j] += (int64_t)params->q[i];
    }
  }

  correct = check_ntt_mod_switch(ct2->ab.b, d0, params);

  if (correct) {
    printf("---------- D0 CORRECT ----------\n");
  } else {
    printf("---------- D0 WRONG ----------\n");
  }

  rns_free(&d0);
}

void check_d1(const ciphertext_t* ct2, const ciphertext_t* ct,
              const rns_params_t* params) {
  printf("---------- VERIFY D1 ----------\n");
  uint8_t correct = 1;

  rns_t d1;
  rns_new(&d1, params);
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l; i++) {
    // D1 = 2A1B1
    for (uint32_t j = 0; j < n; j++) {
      d1[n * i + j] =
          (int64_t)(((__int128_t)ct->ab.b[n * i + j] * ct->ab.a[n * i + j]) %
                    (int64_t)params->q[i]);
      d1[n * i + j] =
          (int64_t)(((__int128_t)2 * d1[n * i + j]) % (int64_t)params->q[i]);
    }
    poly_intt(&d1[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&d1[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    for (uint32_t j = 0; j < n; j++) {
      if (d1[n * i + j] < 0) d1[n * i + j] += (int64_t)params->q[i];
    }
  }

  correct = check_ntt_mod_switch(ct2->ab.a, d1, params);

  if (correct) {
    printf("---------- D1 CORRECT ----------\n");
  } else {
    printf("---------- D1 WRONG ----------\n");
  }

  rns_free(&d1);
}

void check_conjugation(const vec_t* v1, const vec_t* v2, uint32_t m) {
  printf("---------- VERIFY CONJUGATION ----------\n");
  uint8_t correct = 1;

  for (uint32_t i = 0; i < m; i++) {
    // TODO: Better metric
    double diff = v1->rel[i] - v2->rel[i];
    if ((diff > 0.5) | (diff < -0.5)) {
      correct = 0;
      printf("v1.rel[%d]: %lf, v2.rel[%d]: %lf\n", i, v1->rel[i], i,
             v2->rel[i]);
    }
    diff = v1->img[i] + v2->img[i];
    // TODO: Better metric
    if ((diff > 0.5) | (diff < -0.5)) {
      correct = 0;
      printf("v1.img[%d]: %lf, v2.img[%d]: %lf\n", i, -v1->img[i], i,
             v2->img[i]);
    }
  }

  if (correct) {
    printf("---------- CONJUGATION CORRECT ----------\n");
  } else {
    printf("---------- CONJUGATION WRONG ----------\n");
  }
}

void check_inner_product(const vec_t* v1, const vec_t* v2, uint32_t m,
                         uint32_t n) {
  printf("---------- VERIFY INNER PRODUCT ----------\n");
  uint8_t correct = 1;

  uint32_t ratio = n / (2 * m);
  double inner_rel1 = 0, inner_img1 = 0;
  double inner_rel2 = v2->rel[0] / ratio, inner_img2 = v2->img[0] / ratio;

  for (uint32_t i = 0; i < m; i++) {
    inner_rel1 += rel_mul(v1->rel[i], v1->img[i], v1->rel[i], v1->img[i]);
    inner_img1 += img_mul(v1->rel[i], v1->img[i], v1->rel[i], v1->img[i]);
  }

  // TODO: Better metric
  double diff = inner_rel1 - inner_rel2;
  if ((diff > 0.5) | (diff < -0.5)) {
    correct = 0;
    printf("inner1_rel: %lf, inner2_rel: %lf\n", inner_rel1, inner_rel2);
  }
  diff = inner_img1 - inner_img2;
  // TODO: Better metric
  if ((diff > 0.5) | (diff < -0.5)) {
    correct = 0;
    printf("inner1_img: %lf, inner2_img: %lf\n", inner_img1, inner_img2);
  }

  if (correct) {
    printf("---------- INNER PRODUCT CORRECT ----------\n");
  } else {
    printf("---------- INNER PRODUCT WRONG ----------\n");
  }
}

#define MIN(A, B) (((A) > (B)) ? (B) : (A))

void check_linear_transform(const vec_t* v1, const vec_t* v2, const vec_t* mat,
                            uint32_t m) {
  printf("---------- VERIFY LINEAR TRANSFORM ----------\n");
  uint8_t correct = 1;
  vec_t v3;

  vec_new(&v3, m);

  for (uint32_t i = 0; i < m; i++) {
    v3.rel[i] = 0;
    v3.img[i] = 0;
  }

  for (uint32_t d = 0; d < m; d++) {
    for (uint32_t i = 0; i < m; i++) {
      uint32_t j = (i + d) % m;
      v3.rel[i] +=
          rel_mul(v1->rel[i], v1->img[i], mat[i].rel[j], mat[i].img[j]);
      v3.img[i] +=
          img_mul(v1->rel[i], v1->img[i], mat[i].rel[j], mat[i].img[j]);
    }
  }

  // TODO: Better metric
  for (uint32_t i = 0; i < m; i++) {
    double diff = v3.rel[i] - v2->rel[i];
    if ((diff > 1.0) | (diff < -1.0)) {
      correct = 0;
      printf("v_rel[%d]: %lf, v2_rel[%d]: %lf\n", i, v3.rel[i], i, v2->rel[i]);
    }
    diff = v3.img[i] - v2->img[i];
    if ((diff > 1.0) | (diff < -1.0)) {
      correct = 0;
      printf("v_img[%d]: %lf, v2_img[%d]: %lf\n", i, v3.img[i], i, v2->img[i]);
    }
  }

  vec_free(&v3);
  if (correct) {
    printf("---------- LINEAR TRANSFORM CORRECT ----------\n");
  } else {
    printf("---------- LINEAR TRANSFORM WRONG ----------\n");
  }
}

void check_hybrid_multiplication(xof_state_t* xof_state,
                                 const rns_params_t* params) {
  printf("---------- VERIFY HYBRID MULTIPLICATION ----------\n");
  rns_t d, d2, d3;
  uint32_t n = params->n;
  uint8_t correct = 1;

  rns_new(&d, params);
  rns_new(&d2, params);
  rns_new(&d3, params);

  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(&d[n * i], xof_state, 0, params->q[i], params->R_q[i],
                         n);
    for (uint32_t j = 0; j < n; j++) {
      d[n * i + j] = 1;
    }
  }

  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(&d2[n * i], &d[n * i], n * sizeof(int64_t));
    poly_ntt(&d[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    for (uint32_t j = 0; j < params->k; j++) {
      poly_mulx(&d[n * i], (int64_t)params->q[params->l + j], params->q[i],
                params->qinvmR[i], n);
      poly_mulx(&d[n * i], params->R2mq[i], params->q[i], params->qinvmR[i], n);
    }
  }
  rns_mod_raise(d, params->l, params);
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    poly_ntt(&d[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
  rns_mod_switch(d, d3, params->l, params->l, params->k, params->zRbq,
                 params->yRbq, params);

  for (uint32_t i = 0; i < params->l; i++) {
    poly_intt(&d[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&d[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    for (uint32_t j = 0; j < n; j++) {
      if (d2[n * i + j] != d[n * i + j]) {
        correct = 0;
        printf("d[%d]: %ld, d2[%d]: %ld\n", n * i + j, d[n * i + j], n * i + j,
               d2[n * i + j]);
      }
    }
  }

  rns_free(&d);
  rns_free(&d2);
  rns_free(&d3);
  if (correct) {
    printf("---------- HYBRID MULTIPLICATION CORRECT ----------\n");
  } else {
    printf("---------- HYBRID MULTIPLICATION WRONG ----------\n");
  }
}

int main() {
  // Parameters Construction
  uint32_t l = 1;
  uint32_t n2 = 16;
  uint32_t q0_bits = 60;
  uint32_t scale_bits = 40;

  ckks_params_t params;

#ifdef __RUN_MODE__
  uint64_t start, end;
#endif

#ifdef __DEBUG_MODE__
  printf("---------- KEY GENERATION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  // 10x less tests because it takes too long
  for (uint32_t i = 0; i < (N_TESTS + 9) / 10; i++) {
    ckks_keygen(&params, q0_bits, scale_bits, n2, l, 0, 0, NULL, 0);
    ckks_relinearization_keygen(&params);
    ckks_conjugation_keygen(&params);
    ckks_public_keygen(&params);
    ckks_inner_product_keygen(&params);
    ckks_linear_transform_keygen(&params);
  }
#ifdef __DEBUG_MODE__
  check_q_aux(params.rns.q, params.rns.Rmq, params.rns.R2mq, params.rns.R_q,
              params.rns.qinvmR, params.rns._nmq, params.rns.R_nmq,
              params.rns.R2_nmq, params.rns.Rqinvmq, params.rns.zRbq,
              params.rns.yRbq, params.rns.zRqb, params.rns.yRqb, params.rns.l,
              params.rns.k);
  check_roots(params.rns.w, params.rns.q, params.rns.l + params.rns.k,
              params.rns.n);
  check_roots(params.rns.w_1, params.rns.q, params.rns.l + params.rns.k,
              params.rns.n);
  check_keys(&params);
  check_hybrid_multiplication(&params.xof_state, &params.rns);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Key Generation: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Encode
#ifdef __DEBUG_MODE__
  printf("---------- ENCODE ----------\n");
#endif
  vec_t v;
  vec_new(&v, n2);
  for (uint32_t i = 0; i < n2; i++) {
    v.rel[i] = (double)0.25 * i;
    v.img[i] = 0;
  }
  // printf("---------- ENCODE VECTOR ----------\n");
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v_rel[%d]: %lf\n", i, v.rel[i]);
  //   printf("v_img[%d]: %lf\n", i, v.img[i]);
  // }
  poly_t m;
  poly_new(&m, params.rns.n);
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_encode(m, &v, &params);
  }
#ifdef __DEBUG_MODE__
  check_encoded_polynomial(m, (int64_t*)params.rns.q, params.rns.n);
#endif
  // printf("---------- ENCODE POLYNOMIAL ----------\n");
  // for (uint32_t i = 0; i < n; i++) {
  //   printf("m[%d]: %ld\n", i, m[i]);
  // }

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Encode: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

// Encrypt
#ifdef __DEBUG_MODE__
  printf("---------- ENCRYPT ----------\n");
#endif
  ciphertext_t ct;
  ciphertext_new(&ct, &params.rns);
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_encrypt(&ct, m, &params);
  }
#ifdef __DEBUG_MODE__
  check_ciphertext(&ct, &params.rns);
#endif

  // printf("---------- ENCRYPTED B ----------\n");
  // for (uint32_t i = 0; i < l; i++) {
  //   for (uint32_t j = 0; j < params.rns.n; j++) {
  //     printf("b[%d][%d]: %ld\n", i, j, ct.ab.b[params.rns.n * i + j]);
  //   }
  // }

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Encrypt: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

#ifdef __DEBUG_MODE__
  printf("---------- DECRYPT ----------\n");
#endif
  poly_t m2;
  poly_new(&m2, params.rns.n);
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_decrypt(m2, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  check_decoded_polynomial(m, m2, params.rns.scales[ct.level - 1],
                           params.rns.n);
#endif
  // printf("---------- DECRYPTED POLYNOMIAL ----------\n");
  // for (uint32_t i = 0; i < params.rns.n; i++) {
  //   printf("m2[%d]: %ld\n", i, m2[i]);
  // }

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Decrypt: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

#ifdef __DEBUG_MODE__
  printf("---------- DECODE ----------\n");
#endif
  vec_t v2;
  vec_new(&v2, n2);
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_decode(&v2, m2, ct.level, &params);
  }
#ifdef __DEBUG_MODE__
  check_vectors(&v, &v2, n2);
#endif
  // printf("---------- DECODED VECTOR ----------\n");
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
  //   printf("v2_img[%d]: %lf\n", i, v2.img[i]);
  // }

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Decode: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  ciphertext_t ct2;
  vec_t v3;
  vec_new(&v3, n2);
  ciphertext_new(&ct2, &params.rns);
  // Addition
#ifdef __DEBUG_MODE__
  printf("---------- CIPHER-CIPHER ADDITION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_ciphertext_addition(&ct2, &ct, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  for (uint32_t i = 0; i < n2; i++) {
    v3.rel[i] = 2 * v.rel[i];
    v3.img[i] = 2 * v.img[i];
  }
  check_vectors(&v3, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Cipher-Cipher Add: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Multiplication
#ifdef __DEBUG_MODE__
  printf("---------- CIPHER-CIPHER MULTIPLICATION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_ciphertext_multiplication(&ct2, &ct, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  for (uint32_t i = 0; i < n2; i++) {
    v3.rel[i] = rel_mul(v.rel[i], v.img[i], v.rel[i], v.img[i]);
    v3.img[i] = img_mul(v.rel[i], v.img[i], v.rel[i], v.img[i]);
  }
  check_vectors(&v3, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Cipher-Cipher Mul: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Addition
#ifdef __DEBUG_MODE__
  printf("---------- CIPHER-PLAIN ADDITION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_plaintext_addition(&ct2, &ct, m, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  for (uint32_t i = 0; i < n2; i++) {
    v3.rel[i] = 2 * v.rel[i];
    v3.img[i] = 2 * v.img[i];
  }
  check_vectors(&v3, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Cipher-Plain Add: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Multiplication
#ifdef __DEBUG_MODE__
  printf("---------- CIPHER-PLAIN MULTIPLICATION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_plaintext_multiplication(&ct2, &ct, m, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  for (uint32_t i = 0; i < n2; i++) {
    v3.rel[i] = rel_mul(v.rel[i], v.img[i], v.rel[i], v.img[i]);
    v3.img[i] = img_mul(v.rel[i], v.img[i], v.rel[i], v.img[i]);
  }
  check_vectors(&v3, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Cipher-Plain Mul: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Conjugation
#ifdef __DEBUG_MODE__
  printf("---------- CONJUGATION ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_conjugation(&ct2, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  check_conjugation(&v, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Conjugation: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Public Key
#ifdef __DEBUG_MODE__
  printf("---------- PUBLIC KEY ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_server_encrypt(&ct2, m, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  check_vectors(&v, &v2, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Public Key: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Inner-Product
#ifdef __DEBUG_MODE__
  printf("---------- INNER PRODUCT ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_ciphertext_inner_product(&ct2, &ct, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  check_inner_product(&v, &v2, n2, params.rns.n);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Inner Product: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Linear Transform
#ifdef __DEBUG_MODE__
  printf("---------- LINEAR TRANSFORM ----------\n");
#endif
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  vec_t* mat1;
  mat1 = malloc(n2 * sizeof(vec_t));
  for (uint32_t i = 0; i < n2; i++) {
    vec_new(&mat1[i], n2);
    for (uint32_t j = 0; j < n2; j++) {
      mat1[i].rel[j] = 0.000001 * (n2 * i + j);
      mat1[i].img[j] = -0.000001 * (n2 * i + j);
    }
  }
  poly_t* mat2;
  mat2 = malloc(params.rns.n * sizeof(poly_t));
  for (uint32_t i = 0; i < params.rns.n; i++) {
    poly_new(&mat2[i], params.rns.n);
  }
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_linear_transform_encode(mat2, mat1, &params);
    ckks_ciphertext_linear_transform(&ct2, &ct, mat2, &params);
  }
#ifdef __DEBUG_MODE__
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  check_linear_transform(&v, &v2, mat1, n2);
#endif

#ifdef __RUN_MODE__
  end = time_ns();
  printf("Linear Transform: %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Free
  ckks_free(&params);
  ciphertext_free(&ct);
  ciphertext_free(&ct2);
  vec_free(&v);
  vec_free(&v2);
  vec_free(&v3);
  poly_free(&m);
  poly_free(&m2);
  for (uint32_t i = 0; i < n2; i++) {
    vec_free(&mat1[i]);
  }
  free(mat1);
  for (uint32_t i = 0; i < params.rns.n; i++) {
    poly_free(&mat2[i]);
  }
  free(mat2);

#ifdef __DEBUG_MODE__
  printf("---------- END OF TESTS ----------\n");
#endif

  return 0;
}
