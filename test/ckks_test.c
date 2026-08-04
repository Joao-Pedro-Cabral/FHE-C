
#include "ckks.h"

#include <stdio.h>

int main() {
  uint32_t n = 4;

  // Generate KeyGen
  vec w;
  vec_new(&w, n, 1);

  ckks_keygen(&w, n);

  // Encode
  poly m;
  poly_new(&m, n);
  uint32_t int_part = 4;
  vec v;
  vec_new(&v, n / 2, int_part);
  v.rel[0] = (int64_t)(1.1 * (1UL << (63 - int_part)));
  v.img[0] = (int64_t)(4.3 * (1UL << (63 - int_part)));
  v.rel[1] = (int64_t)(3.5 * (1UL << (63 - int_part)));
  v.img[1] = (int64_t)(-1.4 * (1UL << (63 - int_part)));
  v.int_part = int_part;
  int64_t scale = 1024;
  ckks_encode(m, &v, &w, scale, n);

  // Decode
  vec v2;
  vec_new(&v2, n / 2, int_part);
  ckks_decode(&v2, m, &w, scale, n);

  // Debug
  printf("---------- DEBUG ----------\n");
  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v_rel[%d]: %lf\n", i, (double)v.rel[i] / (1UL << (63 - int_part)));
    printf("v_img[%d]: %lf\n", i, (double)v.img[i] / (1UL << (63 - int_part)));
  }

  for (uint32_t i = 0; i < n; i++) {
    printf("m[%d]: %ld\n", i, m[i]);
  }

  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v2_rel[%d]: %ld\n", i, v2.rel[i]);
    printf("v2_img[%d]: %ld\n", i, v2.img[i]);
  }

  // Free
  vec_free(&w);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);

  return 0;
}
