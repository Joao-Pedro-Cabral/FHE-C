
#ifndef __RNS_H__
#define __RNS_H__

#include "poly.h"

typedef poly* rns;  // L polynomials (one per q_i)

void rns_new(rns* a, uint32_t n, uint32_t l);

void rns_free(rns* a, uint32_t l);

void rns_gen_secret_key(rns s, const rns roots, const uint64_t* q, uint32_t n,
                        uint32_t l);

void rns_gen_roots(rns roots, const uint64_t* q, uint32_t n, uint32_t l);

void rns_gen_inverse_roots(rns roots, const uint64_t* q, uint32_t n,
                           uint32_t l);

#endif  // __RNS_H__
