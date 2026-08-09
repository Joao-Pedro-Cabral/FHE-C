#ifndef __XOF_H__
#define __XOF_H__

#include "blake3.h"
#include <stddef.h>
#include <stdint.h>

typedef blake3_hasher xof_state;

#define XOF_RATE BLAKE3_BLOCK_LEN

/* One-stop XOF call */
void xof(uint8_t *output, size_t outlen, const uint8_t *input,
         size_t inlen);

/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void xof_absorb(xof_state *state, const uint8_t *input, size_t inlen);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void xof_squeeze(uint8_t *output, size_t outlen, xof_state *state);

#endif // __XOF_H__
