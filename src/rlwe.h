
#ifndef __RLWE_H__
#define __RLWE_H__

#include "poly.h"
#include "set.h"

typedef poly* rlwe;  // L polynomials (one per q_i)

void rlwe_new(rlwe* a, uint32_t n, uint32_t l);

void rlwe_free(rlwe* a, uint32_t l);

void rlwe_gen_secret(rlwe s, const set eta, const uint64_t* q, uint32_t n,
                     uint32_t l);

#endif  // __RLWE_H__
