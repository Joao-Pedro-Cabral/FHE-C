
#ifndef __RLWE_H__
#define __RLWE_H__

#include "rns.h"

typedef struct {
  rns_t b;  // ciphertext (b = as + deltam + e)
  rns_t a;  // ciphertext (a = -a)
} rlwe_t;

typedef rlwe_t* switch_t;  // L RLWS's for Key Switch (one per q_i)

// RLWE Keys
typedef struct {
  poly_t s_hat;  // secrt key in normal domain
  rns_t s;       // secret key in l NTT domains
  switch_t r;    // relinearization key (l RLWE's)
  switch_t* h;   // log2_n - 1 rotation keys by 5^h, h = 2^i
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
void rlwe_free(rlwe_t* ab, const rns_params_t* params);
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
  @param params: RNS params
*/
void rlwe_gen_secret_key(keys_t* keys, const rns_params_t* params);

/*
  @param keys: keys
  @param params: RNS params
*/
void rlwe_gen_relinearization_key(keys_t* keys, const rns_params_t* params);

/*
  @param keys: keys
  @param params: RNS params
*/
void rlwe_gen_rotation_key(keys_t* keys, const rns_params_t* params);

/*
  @param m: input polynomial
  @param keys: RLWE context
  @param params: RNS params
  @param ab: output encrypted RLWE
*/
void rlwe_plaintext_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
                            const rns_params_t* params);

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
  @param keys: RLWE context
  @param key_switch: output key switch
  @param e: error polynomial (allocated by caller)
  @param m: output polynomial
*/
void rlwe_key_switch_encrypt(switch_t key_switch, const poly_t m, poly_t me,
                             poly_t e, const keys_t* keys,
                             const rns_params_t* params);

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
  @param tmp: temporary
  @param tmp2: temporary 2
  @param ab: output RLWE
*/
void rlwe_key_switch(rlwe_t* ab, const switch_t r, poly_t tmp, poly_t tmp2,
                     uint32_t level, const rns_params_t* params);

#endif  // __RLWE_H__
