
#include "ckks.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "complex.h"
#include "random.h"
#include "utils.h"

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

void check_q_aux(const uint64_t* q, const int64_t* qR, const int64_t* qR2,
                 const int64_t* qRq, const int64_t* qinv, const int64_t* qn,
                 const int64_t* qRn, const int64_t* qR2n, uint32_t l) {
  printf("---------- VERIFY Q AUXILIARIES ROOTS ----------\n");
  for (uint32_t i = 0; i < l; i++) {
    printf(
        "qR: %ld, qR2: %ld, qinv: %ld, qRq: %ld, qn: %ld, qRn: %ld , qR2n: %ld "
        "mod "
        "%ld\n",
        qR[i], qR2[i], qinv[i], qRq[i], qn[i], qRn[i], qR2n[i], q[i]);
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

void check_secret_key(const ckks_params_t* params) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly_t s_hat;
  poly_new(&s_hat, params->rns.n);
  for (uint32_t i = 0; i < params->rns.l; i++) {
    memcpy(s_hat, params->keys.s[i], params->rns.n * sizeof(uint64_t));
    poly_intt(s_hat, params->rns.w_1[i], params->rns.q[i], params->rns.qRq[i],
              params->rns.qinv[i], params->rns.n);
    poly_mulx(s_hat, params->rns.qRn[i], params->rns.q[i], params->rns.qinv[i],
              params->rns.n);  // multiply by R^2/n
    for (uint32_t j = 0; j < params->rns.n; j++) {
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
    rlwe_key_switch_decrypt(s2, params->keys.h[i], &params->keys, &params->rns);
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
  check_relinearization_key(params);
  check_rotation_key(params);
}

void check_ciphertext(const ciphertext_t* ct, const rns_params_t* params) {
  printf("---------- VERIFY CYPHERTEXT ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < params->l; i++) {
    for (uint32_t j = 0; j < params->n; j++) {
      if ((ct->ab.a[i][j] > (int64_t)params->q[i] / 2) |
          (ct->ab.a[i][j] < -(int64_t)params->q[i] / 2)) {
        correct = 0;
        printf("a[%d][%d]: %ld mod %ld\n", i, j, ct->ab.a[i][j], params->q[i]);
      }
      if ((ct->ab.b[i][j] > (int64_t)params->q[i] / 2) |
          (ct->ab.b[i][j] < -(int64_t)params->q[i] / 2)) {
        correct = 0;
        printf("b[%d][%d]: %ld mod %ld\n", i, j, ct->ab.b[i][j], params->q[i]);
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

void check_mod_switch(const rns_params_t* params) {
  printf("---------- VERIFY MOD SWITCH ----------\n");
  uint8_t correct = 1;
  rns_t a, a2;
  poly_t tmp;
  rns_new(&a, params);
  rns_new(&a2, params);
  poly_new(&tmp, params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    randombytes((uint8_t*)a[i], params->n * sizeof(int64_t));
    for (uint32_t j = 0; j < params->n; j++) {
      a[i][j] = a[i][j] % (int64_t)params->q[i];
      if (a[i][j] < 0) a[i][j] += (int64_t)params->q[i];
    }
  }
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(a2[i], a[i], params->n * sizeof(int64_t));
  }
  for (uint32_t i = 0; i < params->l; i++) {
    poly_ntt(a[i], params->w[i], params->q[i], params->qRq[i], params->qinv[i],
             params->n);
  }
  rns_mod_switch(a, tmp, params->l - 1, params);
  for (uint32_t i = 0; i < params->l - 1; i++) {
    poly_intt(a[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_mulx(a[i], params->qRn[i], params->q[i], params->qinv[i], params->n);
    for (uint32_t j = 0; j < params->n; j++) {
      int64_t new_a =
          (int64_t)(((__int128_t)params->q[params->l - 1] * a[i][j] +
                     a[params->l - 1][j]) %
                    (int64_t)params->q[i]);
      if (new_a < 0) new_a += (int64_t)params->q[i];
      int64_t diff = new_a - a2[i][j];
      int64_t error_threshold = 2;  // Only removes one moduli
      if ((diff >= error_threshold) || (diff <= -error_threshold)) {
        correct = 0;
        printf("Incorrect new_a: %ld, a[%d][%d]: %ld\n", new_a, i, j, a2[i][j]);
      }
    }
  }
  poly_free(&tmp);
  rns_free(&a, params);
  rns_free(&a2, params);
  if (correct) {
    printf("---------- MOD SWITCH CORRECT ----------\n");
  } else {
    printf("---------- MOD SWITCH WRONG ----------\n");
  }
}

int main() {
  // Parameters Construction
  uint32_t l = 1;
  uint32_t q0_bits = 62;
  uint32_t scale_bits = 50;

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
    ckks_keygen(&params, q0_bits, scale_bits, l);
    ckks_relinearization_keygen(&params);
    ckks_rotation_keygen(&params);
  }
#ifdef __DEBUG_MODE__
  check_q_aux(params.rns.q, params.rns.qR, params.rns.qR2, params.rns.qRq,
              params.rns.qinv, params.rns.qn, params.rns.qRn, params.rns.qR2n,
              params.rns.l);
  check_roots(params.rns.w, params.rns.q, params.rns.l, params.rns.n);
  check_roots(params.rns.w_1, params.rns.q, params.rns.l, params.rns.n);
  check_keys(&params);
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
  vec_new(&v, params.rns.n / 2);
  for (uint32_t i = 0; i < params.rns.n / 2; i++) {
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
    poly_new(&m, params.rns.n);
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
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ciphertext_new(&ct, &params.rns);
    ckks_encrypt(&ct, m, &params);
  }
#ifdef __DEBUG_MODE__
  check_ciphertext(&ct, &params.rns);
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
    poly_new(&m2, params.rns.n);
    ckks_decrypt(m2, &ct, &params);
  }
#ifdef __DEBUG_MODE__
  check_decoded_polynomial(m, m2, params.rns.scales[ct.level - 1],
                           params.rns.n);
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
    vec_new(&v2, params.rns.n / 2);
    ckks_decode(&v2, m2, ct.level, &params);
  }
#ifdef __DEBUG_MODE__
  check_vectors(&v, &v2, params.rns.n / 2);
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

// Multiplication
#ifdef __DEBUG_MODE__
  printf("---------- CIPHER-CIPHER MULTIPLICATION ----------\n");
#endif
  ciphertext_t ct2;
#ifdef __RUN_MODE__
  start = time_ns();
#endif
  for (uint32_t i = 0; i < N_TESTS; i++) {
    ciphertext_new(&ct2, &params.rns);
    ckks_ciphertext_ciphertext_multiplication(&ct2, &ct, &ct, &params);
  }
  ckks_decrypt(m2, &ct2, &params);
  ckks_decode(&v2, m2, ct2.level, &params);
  for (uint32_t i = 0; i < params.rns.n / 2; i++) {
    double t_rel, t_img;
    t_rel = v.rel[i];
    t_img = v.img[i];
    v.rel[i] = rel_mul(t_rel, t_img, t_rel, t_img);
    v.img[i] = img_mul(t_rel, t_img, t_rel, t_img);
  }
#ifdef __DEBUG_MODE__
  check_vectors(&v, &v2, params.rns.n / 2);
  check_mod_switch(&params.rns);
#endif

  // Free
  ckks_free(&params);
  ciphertext_free(&ct, &params.rns);
  ciphertext_free(&ct2, &params.rns);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);
  poly_free(&m2);

#ifdef __DEBUG_MODE__
  printf("---------- END OF TESTS ----------\n");
#endif

  return 0;
}
