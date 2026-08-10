
#include "ckks.h"

#include <stdio.h>
#include <string.h>

void check_q_aux(const uint64_t* q, const int64_t* qR, const int64_t* qinv,
                 const int64_t* qn, uint32_t l) {
  printf("---------- VERIFY Q AUXILIARIES ROOTS ----------\n");
  for (uint32_t i = 0; i < l; i++) {
    printf("qR: %ld, qinv: %ld, qn: %ld mod %ld\n", qR[i], qinv[i], qn[i],
           q[i]);
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

void check_secret_key(const rns s, const uint64_t* q, const int64_t* qR,
                      const int64_t* qinv, const int64_t* qn, const rns roots,
                      uint32_t l, uint32_t n) {
  printf("---------- VERIFY SECRET KEY ----------\n");
  uint8_t correct = 1;
  poly s_hat, ones;
  poly_new(&s_hat, n);
  poly_new(&ones, n);
  for (uint32_t i = 0; i < n; i++) {
    ones[i] = 1;
  }
  for (uint32_t i = 0; i < l; i++) {
    memcpy(s_hat, s[i], n * sizeof(uint64_t));
    poly_intt(s_hat, roots[i], q[i], qR[i], qinv[i], qn[i], n);
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
  if (correct) {
    printf("---------- SECRET KEY CORRECT ----------\n");
  } else {
    printf("---------- SECRET KEY WRONG ----------\n");
  }
}

void check_polynomials(const poly m, const poly m2, uint64_t scale,
                       uint32_t n) {
  printf("---------- VERIFY POLYNOMIALS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < n; i++) {
    int64_t diff = m[i] - m2[i];
    if ((diff > (int64_t)scale / 2) | (diff < -(int64_t)scale / 2)) {
      correct = 0;
      printf("m[%d]: %lx, m2[%d]: %lx\n", i, m[i], i, m2[i]);
    }
  }
  if (correct) {
    printf("---------- POLYNOMIALS CORRECT ----------\n");
  } else {
    printf("---------- POLYNOMIALS WRONG ----------\n");
  }
}

void check_vectors(const vec* v, const vec* v2, uint32_t n) {
  printf("---------- VERIFY VECTORS ----------\n");
  uint8_t correct = 1;
  for (uint32_t i = 0; i < n; i++) {
    double diff = v->rel[i] - v2->rel[i];
    // TODO: Better metric
    if ((diff > 0.01) | (diff < -0.01)) {
      correct = 0;
      printf("v_rel[%d]: %lf, v2_rel[%d]: %lf\n", i, v->rel[i], i, v2->rel[i]);
    }
    diff = v->img[i] - v2->img[i];
    // TODO: Better metric
    if ((diff > 0.01) | (diff < -0.01)) {
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
  uint32_t n = 128;
  uint32_t l = 2;
  uint64_t q[2] = {534529, 747521};
  uint64_t scale = 10000;

  ckks_params params;

  printf("---------- INIT ----------\n");
  ckks_new(&params, n, l, q, scale);
  check_q_aux(params.q, params.qR, params.qinv, params.qn, params.l);
  check_roots(params.roots, params.q, params.l, params.n);
  check_secret_key(params.s, params.q, params.qR, params.qinv, params.qn,
                   params.roots, params.l, params.n);

  // Encode
  printf("---------- ENCODE ----------\n");
  printf("---------- ENCODE VECTOR ----------\n");
  vec v;
  vec_new(&v, n / 2);
  for (uint32_t i = 0; i < n / 2; i++) {
    v.rel[i] = i / 4;
    v.img[i] = i / 4 + 1;
  }
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v_rel[%d]: %lf\n", i, v.rel[i]);
  //   printf("v_img[%d]: %lf\n", i, v.img[i]);
  // }
  poly m;
  poly_new(&m, n);
  ckks_encode(m, &v, &params);
  // printf("---------- ENCODE POLYNOMIAL ----------\n");
  // for (uint32_t i = 0; i < n; i++) {
  //   printf("m[%d]: %ld\n", i, m[i]);
  // }

  // Encrypt
  printf("---------- ENCRYPT ----------\n");
  cyphertext ct;
  cyphertext_new(&ct, &params);
  ckks_encrypt(&ct, m, &params);
  // printf("---------- ENCRYPTED B ----------\n");
  // for (uint32_t i = 0; i < l; i++) {
  //   for (uint32_t j = 0; j < n; j++) {
  //     printf("b[%d][%d]: %ld\n", i, j, ct.b[i][j]);
  //   }
  // }

  printf("---------- DECRYPT ----------\n");
  poly m2;
  poly_new(&m2, n);
  ckks_decrypt(m2, &ct, &params);
  check_polynomials(m, m2, scale, n);
  // printf("---------- DECRYPTED POLYNOMIAL ----------\n");
  // for (uint32_t i = 0; i < n; i++) {
  //   printf("m2[%d]: %ld\n", i, m2[i]);
  // }

  printf("---------- DECODE ----------\n");
  vec v2;
  vec_new(&v2, n / 2);
  ckks_decode(&v2, m2, &params);
  check_vectors(&v, &v2, n / 2);
  // printf("---------- DECODED VECTOR ----------\n");
  // for (uint32_t i = 0; i < n / 2; i++) {
  //   printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
  //   printf("v2_img[%d]: %lf\n", i, v2.img[i]);
  // }

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
