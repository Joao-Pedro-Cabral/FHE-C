#ifndef __XOF_H__
#define __XOF_H__

#include <stddef.h>
#include <stdint.h>

#include "aes128ctr.h"

typedef aes128ctr_ctx xof_state_t;

#define XOF_RATE AES128CTR_BLOCKBYTES
#define SEED_BYTES (128 / 8)

/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void xof_absorb(xof_state_t *state, const uint8_t input[SEED_BYTES]);

/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void xof_squeeze(uint8_t *output, size_t outlen, xof_state_t *state);

/* Squeeze output out of the sponge using OpenMP.
 *
 * Supports being called multiple times
 */
void xof_squeeze_omp(uint8_t *output, size_t outlen, xof_state_t *state);

#endif  // __XOF_H__
