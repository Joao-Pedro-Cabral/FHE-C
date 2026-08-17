
#include "ciphertext.h"

void ciphertext_new(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_new(&ct->ab, params);
  ct->n = params->n;
  ct->l = params->l;
  ct->level = params->l;
}

void ciphertext_free(ciphertext_t* ct, const rns_params_t* params) {
  rlwe_free(&ct->ab, params);
}

void ciphertext_ciphertext_multiplication(ciphertext_t* ct3, const ciphertext_t* ct1, const ciphertext_t* ct2,
                  const switch_t r, const rns_params_t* params) {
  if (ct1->level != ct2->level) {
    printf("Illegal ciphertext multiplication! Both operands must have same CRT!\n");
  }
}
