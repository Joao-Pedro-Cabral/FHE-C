
#include "ckks.h"

#include "random.h"

void ckks_new(ckks_params* params, uint32_t n, uint32_t l, const uint64_t* q,
              uint64_t scale) {
  params->n = n;
  params->l = l;
  params->q = q;
  params->scale = scale;
  // moduli auxiliares
  rns_gen_q_new(&params->qR2n, &params->qR, &params->qRq, &params->qinv, l);
  rns_gen_q_aux(params->qR2n, params->qR, params->qRq, params->qinv, q, n, l);
  // complex roots
  vec_new(&params->w, n);
  vec_gen_roots(&params->w, n);
  // roots for NTT
  rns_new(&params->roots, n, l);
  rns_gen_roots(params->roots, q, params->qR, params->qinv, n, l);
  // s
  rns_new(&params->s, n, l);
  rns_gen_secret_key(params->s, params->roots, q, params->qRq, params->qinv, n,
                     l);
  // roots for INTT
  rns_gen_inverse_roots(params->roots, q, params->qinv, n, l);
}

void ckks_free(ckks_params* params) {
  rns_gen_q_free(&params->qR2n, &params->qR, &params->qRq, &params->qinv);
  vec_free(&params->w);
  rns_free(&params->roots, params->l);
  rns_free(&params->s, params->l);
}

void ckks_encode(poly m, const vec* v, const ckks_params* params) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, params->n);
  vec_hermitian(&v2, v, params->n);
  // Second step: IFFT
  ifft(&v2, &params->w, params->n);
  // Third step: Rescale
  vector_to_poly(m, &v2, params->scale, params->n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly m, const ckks_params* params) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, params->n);
  poly_to_vector(&v2, m, params->scale, params->n);
  // Second step: FFT
  fft(&v2, &params->w, params->n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, params->n);
  // Free
  vec_free(&v2);
}

void ckks_encrypt(cyphertext* ct, const poly m, const ckks_params* params) {
  // DeltaM + E
  poly e;
  poly_new(&e, params->n);
  poly_cbd24_sampler_add(e, m, params->n);
  // B = -(-A)S + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(ct->a[i], params->q[i], params->qRq[i],
                         params->n);  // -A
    poly_pointwise_mul(ct->b[i], ct->a[i], params->s[i], params->q[i],
                       params->qinv[i], params->n);
    poly_intt(ct->b[i], params->roots[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_rsub(ct->b[i], e, params->q[i], params->qR[i], params->qinv[i],
              params->qR2n[i], params->n);
  }
  poly_free(&e);
}

void ckks_decrypt(poly m, const cyphertext* ct, const ckks_params* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  // B + (-A)S = DeltaM + E
  poly_pointwise_mul(m, ct->a[0], params->s[0], params->q[0], params->qinv[0],
                     params->n);
  poly_intt(m, params->roots[0], params->q[0], params->qRq[0], params->qinv[0],
            params->n);
  poly_add(m, ct->b[0], params->q[0], params->qR[0], params->qinv[0],
           params->qR2n[0], params->n);
}

void cyphertext_new(cyphertext* ct, const ckks_params* params) {
  rns_new(&(ct->b), params->n, params->l);
  rns_new(&(ct->a), params->n, params->l);
  ct->n = params->n;
  ct->l = params->l;
  ct->scale = params->scale;
  ct->int_scale = 1;
  ct->noise_scale = 1;
  ct->encoding_type = 4;  // OpenFHE: CKKS_PACKED_ENCODING
}

void cyphertext_free(cyphertext* ct, const ckks_params* params) {
  rns_free(&(ct->b), params->l);
  rns_free(&(ct->a), params->l);
}
