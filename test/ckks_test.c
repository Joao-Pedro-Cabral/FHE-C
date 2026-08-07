
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
  poly s_hat;
  poly_new(&s_hat, n);
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s_hat, s[i], n * sizeof(uint64_t));
    poly_intt(s_hat, eta[i], q[i], n);
    for (uint32_t j = 0; j < n; j++) {
      if (s_hat[j] * s_hat[j] > 1)
        printf("Invalid secret key! s[%d]: %ld mod %ld\n", j, s_hat[j], q[i]);
    }
  }
  poly_free(&s_hat);
}

int main() {
  // Parameters Construction
  uint32_t n = 4;
  uint32_t l = 2;
  uint64_t q[2] = {17, 41};
  int64_t scale = 1024;

  ckks_params params;

  printf("---------- INIT ----------\n");
  ckks_new(&params, n, l, q, scale);
  // check_eta(params.eta, params.q, params.l, params.n);
  check_secret_key(params.s, params.q, params.eta, params.l, params.n);

  // Encode
  poly m;
  poly_new(&m, n);
  uint32_t int_part = 4;
  vec v;
  vec_new(&v, n / 2, int_part);
  v.rel[0] = 1;
  v.img[0] = 4;
  v.rel[1] = 3;
  v.img[1] = -1;
  v.int_part = int_part;
  ckks_encode(m, &v, &params);

  // Decode
  vec v2;
  vec_new(&v2, n / 2, int_part);
  ckks_decode(&v2, m, &params);

  // Debug
  // printf("---------- DEBUG ----------\n");
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v_rel[%d]: %lf\n", i, v.rel[i]);
  //   printf("v_img[%d]: %lf\n", i, v.img[i]);
  // }

  // for (uint32_t i = 0; i < n; i++) {
  //   printf("m[%d]: %ld\n", i, m[i]);
  // }

  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
  //   printf("v2_img[%d]: %lf\n", i, v2.img[i]);
  // }

  // Free
  ckks_free(&params);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);

  return 0;
}
