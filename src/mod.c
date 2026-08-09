
#include "mod.h"

int64_t mod_inverse(uint64_t a, uint64_t q) {
  int64_t t = 0, new_t = 1;
  uint64_t r = q, new_r = a % q;

  while (new_r != 0) {
    uint64_t quotient = r / new_r;

    int64_t tmp = new_t;
    new_t = t - (int64_t)quotient * new_t;
    t = tmp;

    uint64_t tmp2 = new_r;
    new_r = r - quotient * new_r;
    r = tmp2;
  }

  return t;
}
