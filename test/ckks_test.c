
#include "ckks.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef __RUN_MODE__
#define N_TESTS 100
#elif __VTUNE_MODE__
#define N_TESTS 100
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

void check_q_aux(const uint64_t* q, const int64_t* qR, const int64_t* qRq,
                 const int64_t* qinv, const int64_t* qR2n, uint32_t l) {
  printf("---------- VERIFY Q AUXILIARIES ROOTS ----------\n");
  for (uint32_t i = 0; i < l; i++) {
    printf("qR: %ld, qinv: %ld, qRq: %ld, qR2n: %ld mod %ld\n", qR[i], qinv[i],
           qRq[i], qR2n[i], q[i]);
  }
}

void check_roots(const rns_t roots, const uint64_t* q, uint32_t l, uint32_t n) {
  printf("---------- VERIFY RNS ROOTS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      if ((roots[i][j] >= (int64_t)q[i]) | (roots[i][j] <= -(int64_t)q[i])) {
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

void check_secret_key(const rns_t s, const uint64_t* q, const int64_t* qRq,
                      const int64_t* qinv, const int64_t* qR2n,
                      const rns_t roots, uint32_t l, uint32_t n) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly_t s_hat, ones, invn;
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

void check_cyphertext(const ciphertext_t* ct, const int64_t* q) {
  printf("---------- VERIFY CYPHERTEXT ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < ct->l; i++) {
    for (uint32_t j = 0; j < ct->n; j++) {
      if ((ct->ab.a[i][j] > (int64_t)q[i] / 2) |
          (ct->ab.a[i][j] < -(int64_t)q[i] / 2)) {
        correct = 0;
        printf("a[%d][%d]: %ld mod %ld\n", i, j, ct->ab.a[i][j], q[i]);
      }
      if ((ct->ab.b[i][j] > (int64_t)q[i] / 2) |
          (ct->ab.b[i][j] < -(int64_t)q[i] / 2)) {
        correct = 0;
        printf("b[%d][%d]: %ld mod %ld\n", i, j, ct->ab.b[i][j], q[i]);
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

int main() {
  // Parameters Construction
  uint32_t n = 524288;
  uint32_t l = 3;
  uint64_t q[3] = {2305843009255636993UL, 2305843009322745857UL,
                   2305843009696038913UL};
  uint64_t scale = 10000;

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
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_keygen(&params, n, l, q, scale);
    ckks_relinearization_keygen(&params);
  }
#ifdef __DEBUG_MODE__
  check_q_aux(params.rns_t.q, params.rns_t.qR, params.rns_t.qRq,
              params.rns_t.qinv, params.rns_t.qR2n, params.rns_t.l);
  check_roots(params.rns_t.w, params.rns_t.q, params.rns_t.l, params.rns_t.n);
  check_secret_key(params.keys.s, params.rns_t.q, params.rns_t.qRq,
                   params.rns_t.qinv, params.rns_t.qR2n, params.rns_t.w,
                   params.rns_t.l, params.rns_t.n);
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
  poly_t m;
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_new(&m, n);
    ckks_encode(m, &v, &params);
  }
#ifdef __DEBUG_MODE__
  check_encoded_polynomial(m, (int64_t*)params.rns_t.q, n);
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
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    cyphertext_new(&ct, &params);
    ckks_encrypt(&ct, m, &params);
  }
#ifdef __DEBUG_MODE__
  check_cyphertext(&ct, (int64_t*)params.rns_t.q);
#endif
  // printf("---------- ENCRYPTED B ----------\n");
  // for (uint32_t i = 0; i < l; i++) {
  //   for (uint32_t j = 0; j < n; j++) {
  //     printf("b[%d][%d]: %ld\n", i, j, ct->ab.b[i][j]);
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
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_new(&m2, n);
    ckks_decrypt(m2, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  check_decoded_polynomial(m, m2, scale, n);
#endif
  // printf("---------- DECRYPTED POLYNOMIAL ----------\n");
  // for (uint32_t i = 0; i < n; i++) {
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
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    vec_new(&v2, n / 2);
    ckks_decode(&v2, m2, &params);
  }
#ifdef __DEBUG_MODE__
  check_vectors(&v, &v2, n / 2);
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

  // Free
  ckks_free(&params);
  cyphertext_free(&ct, &params);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);
  poly_free(&m2);

#ifdef __DEBUG_MODE__
  printf("---------- END OF TESTS ----------\n");
#endif

  return 0;
}
