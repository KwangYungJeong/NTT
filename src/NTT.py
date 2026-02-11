
# ==========================================
# Number Theoretic Transform (NTT) Implementation
# ==========================================
# Level 1: Core Algorithm
# Optimized with pre-computation and class-based configuration.
# ==========================================

class NTTContext:
    """
    A context for NTT operations with a specific prime and primitive root.
    Provides optimized transforms by pre-computing twiddle factors and bit-reversal maps.
    """
    
    def __init__(self, mod=469762049, g=3):
        """
        Initializes the NTT context.
        
        Args:
            mod (int): The prime modulus (p = c * 2^k + 1).
            g (int): A primitive root modulo p.
        """
        self.mod = mod
        self.g = g
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
    print("-" * 50)
    print("Enhanced NTT (Level 1) - Step-by-Step Demo")
    print(f"Modulus: {_default_ctx.mod}, Root: {_default_ctx.g}")
    print("-" * 50)

    # Goal: (1 + 2x + 3x^2) * (4 + 5x)
    # Result: 4 + (5+8)x + (10+12)x^2 + 15x^3 = 4 + 13x + 22x^2 + 15x^3
    p1 = [1, 2, 3]
    p2 = [4, 5]
    
    print(f"Goal: Multiply A(x) and B(x)")
    print(f"  A(x) = {p1}")
    print(f"  B(x) = {p2}")
    print("-" * 50)

    # --- Step 1: Determine the size of the operation (N) ---
    # The degree of product P(x) = A(x) * B(x) is deg(A) + deg(B).
    # Therefore, the number of coefficients in the result is (len(A) + len(B) - 1).
    target_len = len(p1) + len(p2) - 1
    
    # NTT (like FFT) requires the buffer size 'N' to be:
    # 1. A power of 2 (for the recursive divide-and-conquer butterfly structure).
    # 2. Greater than or equal to 'target_len' (to avoid "cyclic convolution" or aliasing).
    #    If N < target_len, the higher-degree terms will 'wrap around' and corrupt 
    #    lower-degree terms (e.g., x^N would wrap back to x^0).
    n = 1 << (target_len - 1).bit_length()
    
    print(f"[Step 1] Determine N (Transform Size)")
    print(f"  - Target result length: {target_len} Method 1: (Sum of degrees + 1)")
    print(f"  - Target result length: {target_len} Method 2: len(A) + len(B) - 1")
    print(f"  - Calculated N: {n} (Smallest power of 2 >= {target_len})")
    
    fa = p1 + [0] * (n - len(p1))
    fb = p2 + [0] * (n - len(p2))
    print(f"  - Padded A to size {n}: {fa}")
    print(f"  - Padded B to size {n}: {fb}")

    print(f"\n[Step 2] Forward NTT: Moving to Frequency Domain...")
    ntt(fa, False)
    ntt(fb, False)
    # Now fa and fb are "Point-Value" representations
    print(f"  A in freq domain: {fa[:4]}...") 

    print(f"\n[Step 3] Point-wise Multiply: Multiplication in Freq Domain is O(N)")
    for i in range(n):
        fa[i] = (fa[i] * fb[i]) % _default_ctx.mod
    
    print(f"\n[Step 4] Inverse NTT: Returning to Coefficient Domain...")
    ntt(fa, True)
    
    print(f"\n[Step 5] Trimming: Result valid up to target length {target_len}")
    result = fa[:target_len]
    print(f"  Final Coefficients: {result}")
    
    expected = [4, 13, 22, 15]
    if result == expected:
        print("\nVerification SUCCESS! ✅ (Product: 4 + 13x + 22x^2 + 15x^3)")
    else:
        print(f"\nVerification FAILED! ❌ Expected: {expected}")
    print("-" * 50)
