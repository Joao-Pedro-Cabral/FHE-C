
#ifndef __CIPHERTEXT_H__
#define __CIPHERTEXT_H__

#include "rlwe.h"

// OpenFHE ciphertext format
typedef struct {
  rlwe_t ab;       // RLWE of the ciphertext_t
  uint32_t level;  // current level (1 to l)
} ciphertext_t;

/*
  @param params: RNS params
  @param ct: ciphertext_t to be allocated
*/
void ciphertext_new(ciphertext_t* ct, const rns_params_t* params);

/*
  @param params: RNS params
  @param ct: ciphertext_t to be set free
*/
void ciphertext_free(ciphertext_t* ct, const rns_params_t* params);

/*
  @param ct3: output ciphertext_t
  @param ct1: input ciphertext_t
  @param ct2: input ciphertext_t
  @param r: relinearization key switch
  @param params: RNS params
*/
void ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                          const ciphertext_t* ct1,
                                          const ciphertext_t* ct2,
                                          const switch_t r,
                                          const rns_params_t* params);

#endif  // __CIPHERTEXT_H__
