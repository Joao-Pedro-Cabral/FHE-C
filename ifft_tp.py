import numpy as np


def bit_reverse(x, bits):
    y = 0
    for _ in range(bits):
        y = (y << 1) | (x & 1)
        x >>= 1
    return y


def fft_tp(a, output_in_bitrev=True):
    a = np.array(a, dtype=complex)

    n = len(a)
    k = round(np.log2(n))

    print(n, k)
    omega = np.exp(-2j * np.pi / (2*n))

    a_hat = np.copy(a)

    for l in range(k - 1, -1, -1):
        for i in range(0, 2**(k - 1 - l)):
            # This is for cyclic FFT
            phi = omega ** -bit_reverse(2 ** (k-1-l) + i, k)
            for j in range(i * 2**(l + 1), i * 2**(l + 1) + 2**l):
                t0 = a_hat[j]
                t1 = phi * a_hat[j + 2**l]
                a_hat[j] = t0 + t1
                a_hat[j + 2**l] = t0 - t1

    if not output_in_bitrev:
        return np.array([a_hat[bit_reverse(i, k)] for i in range(n)])
    return a_hat


def ifft_tp(a_hat, input_in_bitrev=True):
    n = len(a_hat)
    k = round(np.log2(n))

    if not input_in_bitrev:
        a_hat = np.array([a_hat[bit_reverse(i, k)] for i in range(n)])
    a_hat = np.array(a_hat, dtype=complex)

    print(n, k)
    omega = np.exp(-2j * np.pi / (2*n))

    a = np.copy(a_hat)

    for l in range(0, k):
        for i in range(0, 2**(k - 1 - l)):
            # This is for cyclic FFT
            phi = omega ** bit_reverse(2 ** (k-1-l) + i, k)
            for j in range(i * 2**(l + 1), i * 2**(l + 1) + 2**l):
                t0 = a[j] + a[j + 2**l]
                t1 = a[j] - a[j + 2**l]
                a[j] = t0
                a[j + 2**l] = phi * t1

    a /= n
    return a


def ifft_tp2(a_hat, input_in_bitrev=True):
    n = len(a_hat)
    k = round(np.log2(n))

    print(n, k)
    omega = np.exp(2j * np.pi / (2*n))

    a = np.copy(a_hat)
    for l in range(k - 1, -1, -1):
        print(f"l: {l}")
        idx = 0
        for i in range(0, 2**(k - 1 - l)):
            for j in range(i * 2**(l + 1), i * 2**(l + 1) + 2**l):
                # Wrong?
                # phi = omega ** (
                #         bit_reverse(5 ** (2 ** l + idx//(2**(k-1-l))), k))
                phi = omega ** bit_reverse(2 ** l + idx//(2**(k-1-l)), k)
                idx += 1
                t0 = a[j] + a[j + 2**l]
                t1 = a[j] - a[j + 2**l]
                a[j] = t0
                a[j + 2**l] = phi * t1

    a /= n
    return a


if __name__ == '__main__':

    x = np.array([1+4j, 3-1j, 3+1j, 1-4j], dtype=complex)

    omega = np.exp(-1j * np.pi / 4)

    print("Omega:")
    print(omega)

    roots = omega ** (2 * np.arange(4) + 1)

    print("Roots:")
    print(roots)

    # Vandermonde matrix
    V = np.vander(roots, N=4, increasing=True).T

    print("V:")
    print(V)

    y = V @ x
    y = y / 4

    print("y")
    print(y)

    y = np.fft.ifft(x)

    print("y fft")
    print(y)

    omega = np.exp(1j * np.pi / 4)

    print("Omega:")
    print(omega)

    roots = omega ** (np.arange(4))

    print("Roots:")
    print(roots)

    y = y * roots

    print("y")
    print(y)

    y = ifft_tp2(x, input_in_bitrev=False)
    #
    print("y_tp")
    print(y)
