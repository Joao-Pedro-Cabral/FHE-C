
#include "ckks.h"

#include <stdio.h>
#include <string.h>

void check_eta(const set eta, const uint64_t* q, uint32_t l, uint32_t n) {
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      printf("eta[%d]: %ld mod %ld\n", j, eta[i][j], q[i]);
    }
  }
}

void check_secret_key(const rlwe s, const uint64_t* q, const set eta,
                      uint32_t l, uint32_t n) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly s_hat;
  poly_new(&s_hat, n);
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s_hat, s[i], n * sizeof(uint64_t));
    poly_intt(s_hat, eta[i], q[i], n);
    for (uint32_t j = 0; j < n; j++) {
      if (s_hat[j] * s_hat[j] > 1) {
        correct = 0;
        printf("Invalid secret key! s[%d]: %ld mod %ld\n", j, s_hat[j], q[i]);
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

int main() {
  // Parameters Construction
  uint32_t n = 4;
  uint32_t l = 2;
  uint64_t q[2] = {337, 353};
  int64_t scale = 40;

  ckks_params params;

  printf("---------- INIT ----------\n");
  ckks_new(&params, n, l, q, scale);
  // check_eta(params.eta, params.q, params.l, params.n);
  check_secret_key(params.s, params.q, params.eta, params.l, params.n);

  // Encode
  printf("---------- ENCODE ----------\n");
  printf("---------- ENCODE VECTOR ----------\n");
  vec v;
  vec_new(&v, n / 2);
  v.rel[0] = 1;
  v.img[0] = 4;
  v.rel[1] = 3;
  v.img[1] = -1;
  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v_rel[%d]: %lf\n", i, v.rel[i]);
    printf("v_img[%d]: %lf\n", i, v.img[i]);
  }
  poly m;
  poly_new(&m, n);
  ckks_encode(m, &v, &params);
  printf("---------- ENCODE POLYNOMIAL ----------\n");
  for (uint32_t i = 0; i < n; i++) {
    printf("m[%d]: %ld\n", i, m[i]);
  }

  // Encrypt
  printf("---------- ENCRYPT ----------\n");
  cyphertext ct;
  cyphertext_new(&ct, &params);
  ckks_encrypt(&ct, m, &params);
  printf("---------- ENCRYPTED B ----------\n");
  for (uint32_t i = 0; i < l; i++) {
    for (uint32_t j = 0; j < n; j++) {
      printf("b[%d][%d]: %ld\n", i, j, ct.b[i][j]);
    }
  }

  printf("---------- DECRYPT ----------\n");
  poly m2;
  poly_new(&m2, n);
  ckks_decrypt(m2, &ct, &params);
  printf("---------- DECRYPTED POLYNOMIAL ----------\n");
  for (uint32_t i = 0; i < n; i++) {
    printf("m2[%d]: %ld\n", i, m[i]);
  }
  for (uint32_t i = 0; i < n; i++) {
    if (m[i] != m2[i])
      printf("m[%d]: %ld differs from m2[%d]: %ld\n", i, m[i], i, m2[i]);
  }

  printf("---------- DECODE ----------\n");
  vec v2;
  vec_new(&v2, n / 2);
  ckks_decode(&v2, m, &params);
  printf("---------- DECODED VECTOR ----------\n");
  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
    printf("v2_img[%d]: %lf\n", i, v2.img[i]);
  }

  // Free
  ckks_free(&params);
  cyphertext_free(&ct, &params);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);
  poly_free(&m2);

  return 0;
}
