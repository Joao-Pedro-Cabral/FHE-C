
#include "xof.h"

void xof(uint8_t* output, size_t outlen, const uint8_t* input, size_t inlen) {
#ifdef __AES_AVX2__
  aes128ctr_prf(output, outlen, input, 0);
#else
  xof_state state;
  blake3_hasher_init(&state);
  blake3_hasher_update(&state, input, inlen);
  blake3_hasher_finalize(&state, output, outlen);
#endif
}

void xof_absorb(xof_state* state, const uint8_t* input, size_t inlen) {
#ifdef __AES_AVX2__
  aes128ctr_init(state, input, 0);
#else
  blake3_hasher_init(state);
  blake3_hasher_update(state, input, inlen);
#endif
}

void xof_squeeze(uint8_t* output, size_t outlen, xof_state* state) {
#ifdef __AES_AVX2__
  aes128ctr_squeezeblocks(output, outlen / XOF_RATE, state);
#else
  blake3_hasher_finalize(state, output, outlen);
#endif
}
