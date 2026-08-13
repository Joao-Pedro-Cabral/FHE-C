
#include "ckks.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef __RUN_MODE__
#define N_TESTS 100
static uint64_t time_ns(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}
#endif

void check_q_aux(const uint64_t* q, const int64_t* qR, const int64_t* qRq,
                 const int64_t* qinv, const int64_t* qR2n, uint32_t l) {
  printf("---------- VERIFY Q AUXILIARIES ROOTS ----------\n");
  for (uint32_t i = 0; i < l; i++) {
    printf("qR: %ld, qinv: %ld, qRq: %ld, qR2n: %ld mod %ld\n", qR[i], qinv[i],
           qRq[i], qR2n[i], q[i]);
  }
}

void check_roots(const rns roots, const uint64_t* q, uint32_t l, uint32_t n) {
  printf("---------- VERIFY RNS ROOTS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      if ((roots[i][j] > (int64_t)q[i] / 2) |
          (roots[i][j] < -(int64_t)q[i] / 2)) {
        correct = 0;
        printf("roots[%d]: %ld mod %ld\n", j, roots[i][j], q[i]);
      }
    }
  }
  if (correct) {
    printf("---------- RNS ROOTS CORRECT ----------\n");
  } else {
    printf("---------- RNS ROOTS WRONG ----------\n");
  }
}

void check_secret_key(const rns s, const uint64_t* q, const int64_t* qRq,
                      const int64_t* qinv, const int64_t* qR2n, const rns roots,
                      uint32_t l, uint32_t n) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly s_hat, ones, invn;
  poly_new(&s_hat, n);
  poly_new(&ones, n);
  poly_new(&invn, n);
  for (uint32_t i = 0; i < n; i++) {
    ones[i] = 1;
  }
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s_hat, s[i], n * sizeof(uint64_t));
    poly_intt(s_hat, roots[i], q[i], qRq[i], qinv[i], n);
    for (uint32_t j = 0; j < n; j++) {
      invn[j] = qR2n[i];
    }
    poly_pointwise_mul(s_hat, s_hat, invn, q[i], qinv[i],
                       n);  // multiply by R^2/n
    poly_pointwise_mul(s_hat, s_hat, ones, q[i], qinv[i],
                       n);  // remove extra montgomery factor
    for (uint32_t j = 0; j < n; j++) {
      if (s_hat[j] * s_hat[j] > 1) {
        correct = 0;
        printf("Invalid secret key! s[%d]: %ld mod %ld\n", j, s_hat[j], q[i]);
      }
    }
  }
  poly_free(&s_hat);
  poly_free(&ones);
  poly_free(&invn);
  if (correct) {
    printf("---------- SECRET KEY CORRECT ----------\n");
  } else {
    printf("---------- SECRET KEY WRONG ----------\n");
  }
}

void check_cyphertext(const cyphertext* ct, const int64_t* q) {
  printf("---------- VERIFY CYPHERTEXT ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < ct->l; i++) {
    for (uint32_t j = 0; j < ct->n; j++) {
      if ((ct->a[i][j] > (int64_t)q[i] / 2) |
          (ct->a[i][j] < -(int64_t)q[i] / 2)) {
        correct = 0;
        printf("a[%d][%d]: %ld mod %ld\n", i, j, ct->a[i][j], q[i]);
      }
      if ((ct->b[i][j] > (int64_t)q[i] / 2) |
          (ct->b[i][j] < -(int64_t)q[i] / 2)) {
        correct = 0;
        printf("b[%d][%d]: %ld mod %ld\n", i, j, ct->b[i][j], q[i]);
      }
    }
  }
  if (correct) {
    printf("---------- CYPHERTEXT CORRECT ----------\n");
  } else {
    printf("---------- CYPHERTEXT WRONG ----------\n");
  }
}

void check_encoded_polynomial(const poly m, const int64_t* q, uint32_t n) {
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

void check_decoded_polynomial(const poly m, const poly m2, uint64_t scale,
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

void check_vectors(const vec* v, const vec* v2, uint32_t n) {
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

int main() {
  // Parameters Construction
  uint32_t n = 524288;
  uint32_t l = 3;
  uint64_t q[3] = {2305843009255636993UL, 2305843009322745857UL,
                   2305843009696038913};
  uint64_t scale = 10000;

  ckks_params params;

#ifdef __RUN_MODE__
  uint64_t start, end;
#endif

  printf("---------- KEY GENERATION ----------\n");
#ifdef __RUN_MODE__
  start = time_ns();
  for (uint32_t i = 0; i < N_TESTS; i++) {
#endif
    ckks_keygen(&params, n, l, q, scale);
#ifndef __RUN_MODE__
    check_q_aux(params.q, params.qR, params.qRq, params.qinv, params.qR2n,
                params.l);
    check_roots(params.roots, params.q, params.l, params.n);
    check_secret_key(params.s, params.q, params.qRq, params.qinv, params.qR2n,
                     params.roots, params.l, params.n);
#endif

#ifdef __RUN_MODE__
  }
  end = time_ns();
  printf("elapsed time = %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Encode
  printf("---------- ENCODE ----------\n");
  vec v;
  vec_new(&v, n / 2);
  for (uint32_t i = 0; i < n / 2; i++) {
    v.rel[i] = i / 4;
    v.img[i] = i / 4 + 1;
  }
  // printf("---------- ENCODE VECTOR ----------\n");
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v_rel[%d]: %lf\n", i, v.rel[i]);
  //   printf("v_img[%d]: %lf\n", i, v.img[i]);
  // }
  poly m;
#ifdef __RUN_MODE__
  for (uint32_t i = 0; i < N_TESTS; i++) {
    start = time_ns();
#endif
    poly_new(&m, n);
    ckks_encode(m, &v, &params);
#ifndef __RUN_MODE__
    check_encoded_polynomial(m, (int64_t*)params.q, n);
#endif
    // printf("---------- ENCODE POLYNOMIAL ----------\n");
    // for (uint32_t i = 0; i < n; i++) {
    //   printf("m[%d]: %ld\n", i, m[i]);
    // }

#ifdef __RUN_MODE__
  }
  end = time_ns();
  printf("elapsed time = %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Encrypt
  printf("---------- ENCRYPT ----------\n");
  cyphertext ct;
#ifdef __RUN_MODE__
  start = time_ns();
  for (uint32_t i = 0; i < N_TESTS; i++) {
#endif
    cyphertext_new(&ct, &params);
    ckks_encrypt(&ct, m, &params);
#ifndef __RUN_MODE__
    check_cyphertext(&ct, (int64_t*)params.q);
#endif
    // printf("---------- ENCRYPTED B ----------\n");
    // for (uint32_t i = 0; i < l; i++) {
    //   for (uint32_t j = 0; j < n; j++) {
    //     printf("b[%d][%d]: %ld\n", i, j, ct.b[i][j]);
    //   }
    // }

#ifdef __RUN_MODE__
  }
  end = time_ns();
  printf("elapsed time = %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  printf("---------- DECRYPT ----------\n");
  poly m2;
#ifdef __RUN_MODE__
  start = time_ns();
  for (uint32_t i = 0; i < N_TESTS; i++) {
#endif
    poly_new(&m2, n);
    ckks_decrypt(m2, &ct, &params);
#ifndef __RUN_MODE__
    check_decoded_polynomial(m, m2, scale, n);
#endif
    // printf("---------- DECRYPTED POLYNOMIAL ----------\n");
    // for (uint32_t i = 0; i < n; i++) {
    //   printf("m2[%d]: %ld\n", i, m2[i]);
    // }

#ifdef __RUN_MODE__
  }
  end = time_ns();
  printf("elapsed time = %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  printf("---------- DECODE ----------\n");
  vec v2;
#ifdef __RUN_MODE__
  start = time_ns();
  for (uint32_t i = 0; i < N_TESTS; i++) {
#endif
    vec_new(&v2, n / 2);
    ckks_decode(&v2, m2, &params);
#ifndef __RUN_MODE__
    check_vectors(&v, &v2, n / 2);
#endif
    // printf("---------- DECODED VECTOR ----------\n");
    // for (uint32_t i = 0; i < n / 2; i++) {
    //   printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
    //   printf("v2_img[%d]: %lf\n", i, v2.img[i]);
    // }

#ifdef __RUN_MODE__
  }
  end = time_ns();
  printf("elapsed time = %.3f ms\n", (double)(end - start) / 1000.0);
#endif

  // Free
  ckks_free(&params);
  cyphertext_free(&ct, &params);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);
  poly_free(&m2);

  printf("---------- END OF TESTS ----------\n");

  return 0;
}
