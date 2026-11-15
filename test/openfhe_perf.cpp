
#include <chrono>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "openfhe.h"

using namespace lbcrypto;

// --------------------------------------------------
// Parameters
// --------------------------------------------------

#ifndef CKKS_N
#define CKKS_N 512
#endif

#ifndef CKKS_N2
#define CKKS_N2 16
#endif

#ifndef CKKS_L
#define CKKS_L 4
#endif

#ifndef CKKS_N_TESTS
#define CKKS_N_TESTS 1000
#endif

static const uint32_t q0_bits = 60;
static const uint32_t scale_bits = 30;

// --------------------------------------------------
// Timing
// --------------------------------------------------

static uint64_t time_ns() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

static void print_perf(const char* name, uint64_t start, uint64_t end,
                       uint32_t iterations) {
  double total = (double)(end - start) / 1000.0;

  std::cout << name << ") Total: " << std::fixed << std::setprecision(3)
            << total << " us, AVG: " << total / iterations << " us"
            << std::endl;
}

// --------------------------------------------------
// CKKS Key Generation
// --------------------------------------------------

static KeyPair<DCRTPoly> perf_keygen(CryptoContext<DCRTPoly> cc) {
  KeyPair<DCRTPoly> keys;

  uint64_t start = time_ns();

  for (uint32_t i = 0; i < (CKKS_N_TESTS + 9) / 10; i++) {
    keys = cc->KeyGen();
  }

  uint64_t end = time_ns();

  print_perf("KeyGen", start, end, (CKKS_N_TESTS + 9) / 10);

  return keys;
}

// --------------------------------------------------
// Relinearization Key Generation
// --------------------------------------------------

static void perf_relinearization_keygen(CryptoContext<DCRTPoly> cc,
                                        PrivateKey<DCRTPoly> sk) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    cc->EvalMultKeyGen(sk);
  }

  uint64_t end = time_ns();

  print_perf("Relinearization KeyGen", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// Conjugation Key Generation
// --------------------------------------------------

static void perf_conjugation_keygen(CryptoContext<DCRTPoly> cc,
                                    PrivateKey<DCRTPoly> sk) {
  uint64_t start = time_ns();

  uint32_t indexConj = 2 * cc->GetRingDimension() - 1;
  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    cc->EvalAutomorphismKeyGen(sk, {indexConj});
  }

  uint64_t end = time_ns();

  print_perf("Conjugation KeyGen", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// Inner Product / Summation Key Generation
// --------------------------------------------------

static void perf_inner_product_keygen(CryptoContext<DCRTPoly> cc,
                                      PrivateKey<DCRTPoly> sk) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    cc->EvalSumKeyGen(sk);
  }

  uint64_t end = time_ns();

  print_perf("Summation KeyGen", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// CKKS Encoding
// --------------------------------------------------

static void perf_encode(CryptoContext<DCRTPoly> cc,
                        const std::vector<double>& values, Plaintext& pt) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    pt = cc->MakeCKKSPackedPlaintext(values);
  }

  uint64_t end = time_ns();

  print_perf("Encode", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// CKKS Decoding
// --------------------------------------------------

static void perf_decode(Plaintext pt) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    auto values = pt->GetCKKSPackedValue();

    (void)values;
  }

  uint64_t end = time_ns();

  print_perf("Decode", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// CKKS Encryption
// --------------------------------------------------

static void perf_encrypt(CryptoContext<DCRTPoly> cc, PublicKey<DCRTPoly> pk,
                         Plaintext pt, Ciphertext<DCRTPoly>& ct) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    ct = cc->Encrypt(pk, pt);
  }

  uint64_t end = time_ns();

  print_perf("Encrypt", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// CKKS Decryption
// --------------------------------------------------

static void perf_decrypt(CryptoContext<DCRTPoly> cc, PrivateKey<DCRTPoly> sk,
                         Ciphertext<DCRTPoly> ct) {
  uint64_t start = time_ns();

  Plaintext pt;

  for (uint32_t i = 0; i < CKKS_N_TESTS; i++) {
    cc->Decrypt(sk, ct, &pt);
  }

  uint64_t end = time_ns();

  print_perf("Decrypt", start, end, CKKS_N_TESTS);
}

// --------------------------------------------------
// Main
// --------------------------------------------------

int main() {
  // ----------------------------------------------
  // CKKS Parameters
  // ----------------------------------------------

  CCParams<CryptoContextCKKSRNS> parameters;

  parameters.SetMultiplicativeDepth(CKKS_L);
  parameters.SetScalingModSize(scale_bits);
  parameters.SetFirstModSize(q0_bits);

  parameters.SetBatchSize(CKKS_N2);
  parameters.SetRingDim(CKKS_N);

  // Experimental parameters.
  // Not intended for production security.
  parameters.SetSecurityLevel(HEStd_NotSet);

  CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

  // ----------------------------------------------
  // Enable CKKS Features
  // ----------------------------------------------

  cc->Enable(PKE);
  cc->Enable(KEYSWITCH);
  cc->Enable(LEVELEDSHE);
  cc->Enable(ADVANCEDSHE);

  std::cout << "OpenFHE CKKS Benchmark" << std::endl;

  std::cout << "Ring Dimension: " << cc->GetRingDimension() << std::endl;

  std::cout << "Batch Size: " << CKKS_N2 << std::endl;

  std::cout << "Multiplicative Depth: " << CKKS_L << std::endl;

  std::cout << "Number of Tests: " << CKKS_N_TESTS << std::endl;

  std::cout << "--------------------------------------" << std::endl;

  // ----------------------------------------------
  // Key Generation
  // ----------------------------------------------

  auto keys = perf_keygen(cc);

  // ----------------------------------------------
  // Evaluation Key Generation
  // ----------------------------------------------

  perf_relinearization_keygen(cc, keys.secretKey);

  perf_conjugation_keygen(cc, keys.secretKey);

  perf_inner_product_keygen(cc, keys.secretKey);

  // ----------------------------------------------
  // Input Vector
  // ----------------------------------------------

  std::vector<double> values(CKKS_N2);

  for (uint32_t i = 0; i < CKKS_N2; i++) {
    values[i] = (double)(i + 1) / 10.0;
  }

  // ----------------------------------------------
  // Encoding
  // ----------------------------------------------

  Plaintext pt;

  perf_encode(cc, values, pt);

  // ----------------------------------------------
  // Decoding
  // ----------------------------------------------

  perf_decode(pt);

  // ----------------------------------------------
  // Encryption
  // ----------------------------------------------

  Ciphertext<DCRTPoly> ct;

  perf_encrypt(cc, keys.publicKey, pt, ct);

  // ----------------------------------------------
  // Decryption
  // ----------------------------------------------

  perf_decrypt(cc, keys.secretKey, ct);

  // ----------------------------------------------
  // End
  // ----------------------------------------------

  std::cout << "--------------------------------------" << std::endl;

  std::cout << "Benchmark completed." << std::endl;

  return 0;
}
