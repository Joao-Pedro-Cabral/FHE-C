
#ifndef __SET_H__
#define __SET_H__

#include <stdint.h>

typedef int64_t* roots;  // N integer elements (TODO: For now 2N)
typedef roots* set;      // L integer elements

void set_new(set* a, uint32_t n, uint32_t l);

void set_free(set* a, uint32_t l);

void set_gen_eta(set eta, roots eta0, const uint64_t* q, uint32_t n,
                 uint32_t l);

void set_gen_reverse_eta(set eta, const roots eta0, const uint64_t* q,
                         uint32_t n, uint32_t l);

#endif  // __SET_H__
