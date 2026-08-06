
#include "ckks.h"

#include <stdio.h>

int main() {
  // Parameters Construction
  uint32_t n = 4;
  int64_t scale = 1024;
  vec w;
  vec_new(&w, 2 * n, 1);

  ckks_params params;
  params.l = 2;
  params.n = n;
  params.scale = scale;
  params.w = &w;

  ckks_init(&params);

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
  printf("---------- DEBUG ----------\n");
  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v_rel[%d]: %lf\n", i, v.rel[i]);
    printf("v_img[%d]: %lf\n", i, v.img[i]);
  }

  for (uint32_t i = 0; i < n; i++) {
    printf("m[%d]: %ld\n", i, m[i]);
  }

  for (uint32_t i = 0; i < n / 2; i++) {
    printf("v2_rel[%d]: %lf\n", i, v2.rel[i]);
    printf("v2_img[%d]: %lf\n", i, v2.img[i]);
  }

  // Free
  vec_free(&w);
  vec_free(&v);
  vec_free(&v2);
  poly_free(&m);

  return 0;
}
