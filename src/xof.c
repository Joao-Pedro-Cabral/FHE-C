
#include "xof.h"

void xof(uint8_t* output, size_t outlen, const uint8_t* input, size_t inlen) {
  xof_state state;
  blake3_hasher_init(&state);
  blake3_hasher_update(&state, input, inlen);
  blake3_hasher_finalize(&state, output, outlen);
}

void xof_absorb(xof_state* state, const uint8_t* input, size_t inlen) {
  blake3_hasher_init(state);
  blake3_hasher_update(state, input, inlen);
}

void xof_squeeze(uint8_t* output, size_t outlen, xof_state* state) {
  blake3_hasher_finalize(state, output, outlen);
}
