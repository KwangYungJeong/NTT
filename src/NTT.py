
# ==========================================
# Number Theoretic Transform (NTT) Implementation
# ==========================================
# Prime: p = 469762049 (7 * 2^26 + 1)
# Primitive Root: g = 3
# Max N: 2^26 (approx 67 million)
# ==========================================

MOD = 469762049
G = 3

def power(a, b, m):
    """Computes (a^b) % m using modular exponentiation."""
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = (res * a) % m
        a = (a * a) % m
        b //= 2
    return res

def mod_inverse(n, m):
    """Computes modular multiplicative inverse of n modulo m."""
    return power(n, m - 2, m)

def ntt(a, invert):
    """
    Performs NTT (Number Theoretic Transform) or Inverse NTT.
    
    Args:
        a (list): Input polynomial coefficients. Length must be a power of 2.
        invert (bool): If True, performs Inverse NTT.
    
    Returns:
        list: Transformed coefficients (in-place modification).
    """
    n = len(a)
    
    # 1. Bit-reversal permutation (Iterative ordering)
    # Reorders the array so that butterfly operations can be done iteratively.
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]

    # 2. Butterfly Operations (Cooley-Tukey)
    length = 2
    while length <= n:
        # Calculate root of unity for this level (w_len)
        # w_len = g^((MOD-1)/length) % MOD
        # Note: (MOD-1) // length is the rotation step size in the exponent.
        
        step = (MOD - 1) // length
        if invert:
            # For Inverse NTT, we use g^(-step) which is mod_inverse(g^step)
            # Or simply g^(MOD-1 - step) because g^(MOD-1) = 1.
            # But calculating inverse of w_len is safer/standard.
            w_len = power(G, step, MOD)
            w_len = mod_inverse(w_len, MOD)
        else:
            w_len = power(G, step, MOD)
            
        for i in range(0, n, length):
            w = 1
            for j in range(length // 2):
                # Butterfly:
                # u = a[i+j]
                # v = a[i+j+len/2] * w
                # a[i+j] = u + v
                # a[i+j+len/2] = u - v
                
                u = a[i + j]
                v = (a[i + j + length // 2] * w) % MOD
                
                a[i + j] = (u + v) % MOD
                a[i + j + length // 2] = (u - v + MOD) % MOD
                
                w = (w * w_len) % MOD
        length <<= 1

    # 3. For Inverse NTT, scale by n^(-1)
    if invert:
        n_inv = mod_inverse(n, MOD)
        for i in range(n):
            a[i] = (a[i] * n_inv) % MOD
            
    return a

def multiply_polynomials(a, b):
    """
    Multiplies two polynomials a(x) and b(x) using NTT.
    
    Args:
        a (list): Coefficients of polynomial A (low degree first).
        b (list): Coefficients of polynomial B (low degree first).
        
    Returns:
        list: Coefficients of product polynomial C = A * B.
    """
    # 1. Determine size N (power of 2) >= len(a) + len(b) - 1
    # Example: deg(A)=1 (len 2), deg(B)=1 (len 2) -> deg(A*B)=2 (len 3) -> N=4
    n = 1
    while n < len(a) + len(b) - 1:
        n <<= 1
        
    # 2. Pad vectors with zeros to size N
    # Make copies to avoid modifying original lists
    fa = a[:] + [0] * (n - len(a))
    fb = b[:] + [0] * (n - len(b))
    
    # 3. Applying NTT (Forward Transform)
    # fa and fb are now in Frequency Domain (Point-value representation)
    ntt(fa, False)
    ntt(fb, False)
    
    # 4. Point-wise multiplication
    # C(x_i) = A(x_i) * B(x_i)
    for i in range(n):
        fa[i] = (fa[i] * fb[i]) % MOD
        
    # 5. Inverse NTT
    # Transform back to Coefficient Domain
    ntt(fa, True)
    
    # The result 'fa' now contains coefficients of A*B.
    # The degree is len(a) + len(b) - 2, so relevant length is len(a) + len(b) - 1.
    return fa

if __name__ == "__main__":
    print("-" * 50)
    print(f"NTT Polynomial Multiplication")
    print(f"Prime: {MOD}, Root: {G}")
    print("-" * 50)

    # Example 1: (x + 1)(x + 2) = x^2 + 3x + 2
    # Coefficients: [1, 1] * [2, 1] -> [2, 3, 1]
    poly_a = [1, 1]
    poly_b = [2, 1]
    
    print(f"\nTest 1: (1 + x) * (2 + x)")
    print(f"  Poly A: {poly_a}")
    print(f"  Poly B: {poly_b}")
    
    result = multiply_polynomials(poly_a, poly_b)
    # Expected length: 2 + 2 - 1 = 3
    valid_len = len(poly_a) + len(poly_b) - 1
    print(f"  Result (Full N): {result}")
    print(f"  Result (Trimmed): {result[:valid_len]}")
    
    # Example 2: Larger polynomial
    # (1 + 2x + 3x^2) * (4 + 5x)
    # = 4 + 5x + 8x + 10x^2 + 12x^2 + 15x^3
    # = 4 + 13x + 22x^2 + 15x^3
    print(f"\nTest 2: (1 + 2x + 3x^2) * (4 + 5x)")
    poly_c = [1, 2, 3]
    poly_d = [4, 5]
    print(f"  Poly C: {poly_c}")
    print(f"  Poly D: {poly_d}")
    
    res2 = multiply_polynomials(poly_c, poly_d)
    valid_len2 = len(poly_c) + len(poly_d) - 1
    print(f"  Result (Trimmed): {res2[:valid_len2]}")
    
    print("-" * 50)
