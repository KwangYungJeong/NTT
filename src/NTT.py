
# ==========================================
# Number Theoretic Transform (NTT) Implementation
# ==========================================
# Level 1: Core Algorithm
# Optimized with pre-computation and class-based configuration.
# ==========================================

import sys
import argparse

try:
    from prime_search import get_primitive_root, get_2_adic_valuation, get_root_of_unity
except ImportError:
    # Fallback if prime_search is not in path or has issues
    from sympy.ntheory import primitive_root
    def get_primitive_root(p): return int(primitive_root(p))
    def get_2_adic_valuation(p):
        v, p1 = 0, p-1
        while p1 % 2 == 0: v, p1 = v + 1, p1 // 2
        return v

class NTTContext:
    """
    A context for NTT operations with a specific prime and primitive root.
    Provides optimized transforms by pre-computing twiddle factors and bit-reversal maps.
    """
    
    def __init__(self, mod=469762049):
        """
        Initializes the NTT context.
        
        Args:
            mod (int): The prime modulus (p = c * 2^k + 1).
        """
        # Note: mod check is still done here or via get_primitive_root
        self.mod = mod
        self.g = get_primitive_root(mod)
        
        # 2-adic valuation for safety check
        self.max_k = get_2_adic_valuation(mod)
        
        # Pre-computations cache
        # Pre-computations cache
        self.rev = {}          # size N -> bit-reversal list
        self.roots = {}        # size N -> twiddle factors list
        self.inv_roots = {}    # size N -> inverse twiddle factors list

    def _prepare(self, n):
        """
        Pre-computes bit-reversal mapping and twiddle factors for size n.
        
        Args:
            n (int): The size of the transform (must be a power of 2).
        """
        if n in self.rev:
            return

        # 1. Pre-compute Bit-reversal indices
        rev = [0] * n
        h = n.bit_length() - 1
        for i in range(n):
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (h - 1))
        self.rev[n] = rev

        # 2. Pre-compute Twiddle factors (roots of unity)
        # Instead of calculating power(G, step, MOD) every time, 
        # we pre-calculate w^i for each level.
        roots = [1] * n
        inv_roots = [1] * n
        
        # Primitive n-th root of unity w = g^((mod-1)/n) % mod
        w_n = pow(self.g, (self.mod - 1) // n, self.mod)
        w_inv = pow(w_n, self.mod - 2, self.mod)
        
        curr_w = 1
        curr_inv = 1
        for i in range(n):
            roots[i] = curr_w
            inv_roots[i] = curr_inv
            curr_w = (curr_w * w_n) % self.mod
            curr_inv = (curr_inv * w_inv) % self.mod
            
        self.roots[n] = roots
        self.inv_roots[n] = inv_roots

    def transform(self, a, invert=False):
        """
        Performs Forward or Inverse NTT in-place.
        
        Mathematical Background:
        NTT maps coefficients {a_i} to point-values {A_j} where 
        A_j = sum_{i=0}^{n-1} a_i * (w^j)^i (mod p).
        This is exactly like DFT but replaces e^(2πi/n) with w.
        """
        n = len(a)
        self._prepare(n)
        
        rev = self.rev[n]
        for i in range(n):
            if i < rev[i]:
                a[i], a[rev[i]] = a[rev[i]], a[i]
        
        length = 2
        while length <= n:
            half = length // 2
            # Select relevant roots of unity for this length
            # If full n-th root is w, then length-th root is w^(n/length)
            step = n // length
            target_roots = self.inv_roots[n] if invert else self.roots[n]
            
            for i in range(0, n, length):
                for j in range(half):
                    u = a[i + j]
                    v = (a[i + j + half] * target_roots[j * step]) % self.mod
                    
                    a[i + j] = (u + v) % self.mod
                    a[i + j + half] = (u - v + self.mod) % self.mod
            length <<= 1
            
        if invert:
            n_inv = pow(n, self.mod - 2, self.mod)
            for i in range(n):
                a[i] = (a[i] * n_inv) % self.mod
        
        return a

    def multiply(self, a, b):
        """
        Fast polynomial multiplication using NTT.
        Complexity: O(N log N)
        """
        target_len = len(a) + len(b) - 1
        n = 1 << (target_len - 1).bit_length()
        
        fa = a + [0] * (n - len(a))
        fb = b + [0] * (n - len(b))
        
        self.transform(fa, False)
        self.transform(fb, False)
        
        for i in range(n):
            fa[i] = (fa[i] * fb[i]) % self.mod
            
        self.transform(fa, True)
        return fa[:target_len]

# Global instance for ease of use
_default_ctx = NTTContext()

def ntt(a, invert=False):
    return _default_ctx.transform(a, invert)

def multiply_polynomials(a, b):
    return _default_ctx.multiply(a, b)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enhanced NTT Polynomial Multiplication Demo")
    parser.add_argument("--prime", type=int, default=469762049, help="Prime modulus to use (default: 469762049)")
    args = parser.parse_args()

    try:
        ctx = NTTContext(args.prime)
    except Exception as e:
        print(f"Error initializing NTT: {e}")
        sys.exit(1)

    print("-" * 50)
    print("Enhanced NTT (Level 1) - Dynamic Configuration")
    print(f"Modulus: {ctx.mod}, Calculated Root: {ctx.g}")
    print("-" * 50)

    p1 = [1, 2, 3]
    p2 = [4, 5]
    
    print(f"Goal: Multiply A(x) and B(x)")
    print(f"  A(x) = {p1}")
    print(f"  B(x) = {p2}")
    
    target_len = len(p1) + len(p2) - 1
    n = 1 << (target_len - 1).bit_length()
    
    # Verify if the prime is suitable for this N using helper
    if n > (1 << ctx.max_k):
        print(f"⚠️  Warning: Modulus {ctx.mod} only supports N up to 2^{ctx.max_k}={1 << ctx.max_k}.")
        print(f"   Target N={n} exceeds the 2-adic valuation of this prime.")
        sys.exit(1)

    fa = p1 + [0] * (n - len(p1))
    fb = p2 + [0] * (n - len(p2))
    
    ctx.transform(fa, False)
    ctx.transform(fb, False)
    
    for i in range(n):
        fa[i] = (fa[i] * fb[i]) % ctx.mod
    
    ctx.transform(fa, True)
    
    result = fa[:target_len]
    print(f"\nFinal Coefficients: {result}")
    
    expected = [4, 13, 22, 15]
    if result == expected:
        print("\nVerification SUCCESS! ✅")
    else:
        print(f"\nVerification FAILED! ❌ Expected: {expected}")
    print("-" * 50)
