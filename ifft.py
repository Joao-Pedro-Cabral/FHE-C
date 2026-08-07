import numpy as np


def bit_reverse(x, bits):
    y = 0
    for _ in range(bits):
        y = (y << 1) | (x & 1)
        x >>= 1
    return y


def bit_reverse_permutation(x):
    x = np.asarray(x)
    n = len(x)
    m = int(np.log2(n))

    indices = np.arange(n)
    rev = np.zeros(n, dtype=int)

    for i in range(n):
        b = f"{i:0{m}b}"      # Binary string with m bits
        rev[i] = int(b[::-1], 2)

    return x[rev]


def gs_ifft_negacyclic(X):
    N = len(X)

    # ψ = primitive 2N-th root
    psi = np.exp(1j * np.pi / (2 * N))

    # ω = primitive N-th root
    omega = psi**2

    a = np.array(X, dtype=complex)

    # -------------------------
    # Gauss-Sande DIF IFFT
    # -------------------------
    m = N

    while m > 1:

        half = m // 2
        step = N // m

        for k in range(0, N, m):

            for j in range(half):

                w = omega ** (-j * step)

                u = a[k + j]
                v = a[k + j + half]

                a[k + j] = u + v
                a[k + j + half] = (u - v) * w

        m //= 2

    # -------------------------
    # Bit reversal
    # -------------------------
    bits = int(np.log2(N))
    out = np.zeros(N, dtype=complex)

    for i in range(N):
        out[bit_reverse(i, bits)] = a[i]

    # -------------------------
    # Normalize
    # -------------------------
    out /= N

    # -------------------------
    # Untwist
    # -------------------------
    for i in range(N):
        out[i] *= psi**(-i)

    return out


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

# y = gs_ifft_negacyclic(x)
#
# print("y3")
# print(y)
