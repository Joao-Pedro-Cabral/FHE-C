#ifndef __XOF_H__
#define __XOF_H__

#ifdef __AES_MODE__
#include "aes128ctr.h"
#else
#include "blake3.h"
#endif
#include <stddef.h>
#include <stdint.h>

#ifdef __AES_MODE__
typedef aes128ctr_ctx xof_state_t;

#define XOF_RATE AES128CTR_BLOCKBYTES
#else
typedef blake3_hasher xof_state_t;

#define XOF_RATE BLAKE3_BLOCK_LEN
#endif

/* One-stop XOF call */
void xof(uint8_t *output, size_t outlen, const uint8_t *input, size_t inlen);

/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void xof_absorb(xof_state_t *state, const uint8_t *input, size_t inlen);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void xof_squeeze(uint8_t *output, size_t outlen, xof_state_t *state);

#endif  // __XOF_H__
