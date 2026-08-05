
#ifndef __RLWE_H__
#define __RLWE_H__

#include "poly.h"

typedef poly* rlwe;  // L polynomials (one per q_i)

typedef int64_t* set;  // L integer elements

void gen_eta(set eta, uint32_t n, uint32_t q);

#endif  // __RLWE_H__
