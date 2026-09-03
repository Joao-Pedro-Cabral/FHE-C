

import numpy as np

# g(X) -> g(X^5(h)) mod X^n+1


def x_to_x_5h(p, h, q):
    """
    Compute p(X^(5^h)) in Z[X]/(X^n + 1).

    Parameters
    ----------
    p : np.ndarray
        Polynomial coefficients:
        p[i] is the coefficient of X^i.
        Expected length <= n.
    n : int
        Ring parameter, i.e. modulus is X^n + 1.
    h : int
        Integer exponent. Computes X -> X^(5^h).

    Returns
    -------
    np.ndarray
        Coefficients of p(X^(5^h)) mod (X^n + 1),
        represented using degrees 0,...,n-1.
    """
    n = len(p)
    p = np.asarray(p)
    result = np.zeros(n, dtype=p.dtype)

    power = pow(5, h)

    for i, coeff in enumerate(p):
        e = i * power
        r = e % (2 * n)

        if r < n:
            result[r] += coeff
        else:
            result[r - n] += q - coeff

    return result

# Find a primitive root
# Test all cases: w^n = -1 and w^2n = 1 mod q


def find_primitive_root(n, q):
    """
    Find an element w in Z_q with multiplicative order exactly 2n.

    Requires 2n | (q - 1), since Z_q^* has order q-1.
    """
    if (q - 1) % (2 * n) != 0:
        return -1

    for w in range(2, q):
        if pow(w, 2 * n, q) != 1:
            continue

        # Check that no smaller divisor of 2n gives 1.
        primitive = True

        for d in range(1, 2 * n):
            if (2 * n) % d == 0:
                if pow(w, d, q) == 1:
                    primitive = False
                    break

        if primitive:
            return w

    return -1

# Create base [w^5^0, w^5^1, ..., w^5^(n/2-1), w^(-5^0), ..., w^(-5^(n/2-1))]


def create_vector(w, n, q):
    """
    Create the vector

        [w^(5^0), w^(5^1), ..., w^(5^(n/2-1)),
         w^(-5^0), w^(-5^1), ..., w^(-5^(n/2-1))]

    with all operations performed modulo q.

    Parameters
    ----------
    w : int
        Element of Z_q.
    n : int
        Must be even.
    q : int
        Modulus.

    Returns
    -------
    np.ndarray
        Vector in Z_q.
    """
    half = n // 2

    # Positive powers: w^(5^h)
    first_half = [
        pow(w, pow(5, h), q)
        for h in range(half)
    ]

    # Negative powers: w^(-5^h) = (w^(-1))^(5^h)
    w_inv = pow(w, -1, q)

    second_half = [
        pow(w_inv, pow(5, h), q)
        for h in range(half)
    ]

    return np.array(first_half + second_half, dtype=object)


def vandermonde_rows(v, q):
    """
    Vandermonde matrix where each ROW corresponds to
    one element of v and contains its successive powers.
    """
    n = len(v)

    V = np.empty((n, n), dtype=object)

    for i, x in enumerate(v):
        for j in range(n):
            V[i, j] = pow(int(x), j, q)

    return V


def vandermonde_columns(v, q):
    """
    Vandermonde matrix where each COLUMN corresponds to
    one element of v.

    Before constructing the matrix, reverse the first half
    and second half of v independently.
    """
    n = len(v)
    half = n // 2

    # Reverse each half independently
    v_reversed = (
        list(v[:half][::-1]) +
        list(v[half:][::-1])
    )

    V = np.empty((n, n), dtype=object)

    for j, x in enumerate(v_reversed):
        for i in range(n):
            V[i, j] = pow(int(x), i, q)

    return V


def rotate_halves(v, h):
    """
    Rotate the first and second halves of v independently
    to the left by h positions.

    Example:
        v = [1, 2, 3, 4, 5, 6, 7, 8]
        h = 1

        returns [2, 3, 4, 1, 6, 7, 8, 5]
    """
    v = np.asarray(v)

    n = len(v)

    half = n // 2

    first_half = np.roll(v[:half], -h)
    second_half = np.roll(v[half:], -h)

    return np.concatenate([first_half, second_half])


q = 193
n = 32
h = 2
w = find_primitive_root(n, q)

np.set_printoptions(linewidth=10000,
                    precision=0,
                    suppress=True)

print(w)

v = create_vector(w, n, q)

V = vandermonde_rows(v, q)

V2 = vandermonde_columns(v, q)

a = np.arange(1, n + 1)

a5h = x_to_x_5h(a, h, q)

a5h_print = np.asarray(a5h, dtype=np.int64)

print(a5h_print)

x = np.mod(V @ a, q)

x2 = rotate_halves(x, h)

x3 = x2[::-1]

n_inv = pow(n, -1, q)

y = np.mod(np.mod(V2 @ x3, q)*n_inv, q)

y_print = np.asarray(y, dtype=np.int64)

print(y_print)
