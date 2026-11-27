
#include "xof.h"

void xof_absorb(xof_state_t* state, const uint8_t input[SEED_BYTES]) {
  aes128ctr_init(state, input, 0);
}

void xof_squeeze(uint8_t* output, size_t outlen, xof_state_t* state,
                 uint64_t domain_separator) {
  aes128ctr_squeeze(output, outlen, state, domain_separator);
}

void xof_squeeze_omp(uint8_t* output, size_t outlen, xof_state_t* state,
                     uint64_t domain_separator) {
  aes128ctr_squeeze_omp(output, outlen, state, domain_separator);
}

void xof_state_update(xof_state_t* state, uint64_t blocks) {
  aes128ctr_state_update(state, blocks);
}
