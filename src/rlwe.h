
#ifndef __RLWE_H__
#define __RLWE_H__

#include <stdbool.h>

#include "rns.h"

#define N1N2_RATIO 2

typedef struct {
  rns_t b;  // ciphertext (b = as + deltam + e)
  rns_t a;  // ciphertext (a = -a)
} rlwe_t;

typedef rlwe_t* switch_t;  // L RLWS's for Key Switch (one per q_i)

// RLWE Keys
typedef struct {
  poly_t s_hat;  // secrt key in normal domain
  rns_t s;       // secret key in l NTT domains
  switch_t r;    // relinearization key (l+k RLWE's)
  rlwe_t pk;     // public key
  switch_t cg;   // conjugation key
  switch_t* ip;  // log2_n - 1 inner product rotation keys by 5^h, h = 2^i
  switch_t* lt;  // n1+n2 linear transform rotation keys by 5^h
} keys_t;

/*
  @param params: RNS params
  @param ab: rlwe_t to be allocated
*/
void rlwe_new(rlwe_t* ab, const rns_params_t* params);

/*
  @param params: RNS params
  @param ab: rlwe_t to be set free
*/
void rlwe_free(rlwe_t* ab);
/*
  @param params: RNS params
  @param keys: rlwe_t context to be allocated
*/
void rlwe_keys_new(keys_t* keys, const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: rlwe_t context to be set free
*/
void rlwe_keys_free(keys_t* keys, const rns_params_t* params);

/*
  @param keys: keys
  @param xof_state: XOF state
  @param params: RNS params
*/
void rlwe_gen_secret_key(keys_t* keys, xof_state_t* xof_state,
                         const rns_params_t* params);

/*
  @param keys: keys
  @param xof_state: XOF state
  @param params: RNS params
*/
void rlwe_gen_relinearization_key(keys_t* keys, xof_state_t* xof_state,
                                  const rns_params_t* params);

/*
  @param m: input polynomial
  @param keys: RLWE context
  @param params: RNS params
  @param ab: output encrypted RLWE
*/
void rlwe_plaintext_encrypt(rlwe_t* ab, xof_state_t* xof_state, const poly_t m,
                            const keys_t* keys, const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                            const rns_params_t* params);

/*
  @param params: RNS params
  @param xof_state: XOF state
  @param tmp, tmp2: temporary polynomials
  @param keys: RLWE context
  @param key_switch: output key switch
  @param e: error polynomial (allocated by caller)
  @param m: output polynomial
  @param params: RNS params
*/
void rlwe_key_switch_encrypt(switch_t key_switch, xof_state_t* xof_state,
                             const poly_t m, poly_t tmp, poly_t tmp2, poly_t e,
                             const keys_t* keys, const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: RLWE context
  @param key_switch: input key switch
  @param m: output polynomial
*/
void rlwe_key_switch_decrypt(poly_t m, const switch_t key_switch,
                             const keys_t* keys, const rns_params_t* params);

/*
  @param params: RNS params
  @param r: input key switch
  @param level: current level
  @param a: rns A
  @param tmp: temporary
  @param ab: output RLWE
*/
void rlwe_key_switch(rlwe_t* ab, const rns_t a, const switch_t r, poly_t tmp,
                     poly_t tmp2, rns_t tmp_r, uint32_t level,
                     const rns_params_t* params);

/*
  @param params: RNS params
  @param xof_state: XOF state
  @param pk: public key
  @param m: plaintext
  @param ab: output RLWE
*/
void rlwe_plaintext_encrypt_server(rlwe_t* ab, xof_state_t* xof_state,
                                   const poly_t m, const rlwe_t pk,
                                   const rns_params_t* params);

#endif  // __RLWE_H__
