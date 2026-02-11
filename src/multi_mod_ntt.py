
# ==========================================
# Multi-Modulus NTT (CRT-NTT) Implementation
# ==========================================
# Level 1 + Level 2 Integration
# Enables polynomial multiplication with large coefficients
# by combining results from multiple NTT-friendly primes.
# ==========================================

from NTT import NTTContext
from CRT import mrc_crt

class MultiModNTT:
    """
    Coordinates NTT multiplication across multiple moduli
    and reconstructs the results using the Chinese Remainder Theorem.
    """
    
    def __init__(self, moduli_configs):
        """
        Args:
            moduli_configs (list): List of (prime, root) tuples.
        """
        self.contexts = [NTTContext(m, g) for m, g in moduli_configs]
        self.moduli = [m for m, g in moduli_configs]

    def multiply(self, a, b):
        """
        Multiplies polynomials a and b using multiple NTT moduli.
        """
        # 1. Perform NTT multiplication for each modulus
        results_per_mod = []
        for ctx in self.contexts:
            res = ctx.multiply(a, b)
            results_per_mod.append(res)
            
        # 2. Reconstruct each coefficient using CRT
        # results_per_mod is a list of results: [res_m1, res_m2, res_m3]
        # We need to take res_m1[i], res_m2[i], res_m3[i] and apply CRT.
        
        final_result = []
        num_coeffs = len(results_per_mod[0])
        
        for i in range(num_coeffs):
            remainders = [res[i] for res in results_per_mod]
            # Use mrc_crt for reconstruction
            # Note: mrc_crt in CRT.py has lots of prints, 
            # for library use we might want a silent version.
            # For this demo, let's just use it.
            val = mrc_crt(self.moduli, remainders)
            final_result.append(val)
            
        return final_result

if __name__ == "__main__":
    # Selected 3 Primes from prime_search.py
    configs = [
        (2013265921, 31),
        (469762049, 3),
        (1811939329, 13)
    ]
    
    mm_ntt = MultiModNTT(configs)
    
    print("-" * 60)
    print("Multi-Modulus NTT Demo (Handling Large Coefficients)")
    print("-" * 60)
    
    # Case: Resulting coefficients exceed the limit of a single prime (~2^31)
    # Let's use large inputs. (e.g., 2^30 * 2^30 for some terms)
    # p1 = [2^30, 2^30], p2 = [2^30, 2^30]
    # (2^30 + 2^30 x) * (2^30 + 2^30 x) = 2^60 + 2*2^60 x + 2^60 x^2
    large_val = 1 << 30
    p1 = [large_val, large_val]
    p2 = [large_val, large_val]
    
    print(f"Polynomial A: {p1}")
    print(f"Polynomial B: {p2}")
    
    # 1. Single Modulus NTT (Should fail due to overflow/wrap-around)
    single_ctx = NTTContext(469762049, 3) # Rank 2 prime
    single_res = single_ctx.multiply(p1, p2)
    print(f"\n[Single Modulus res (MOD=469762049)]")
    print(f"  Result: {single_res}")
    
    # 2. Multi-Modulus NTT
    print(f"\n[Multi-Modulus res (3 Primes)]")
    multi_res = mm_ntt.multiply(p1, p2)
    print(f"  Result: {multi_res}")
    
    # Verification
    expected = [large_val * large_val, 2 * large_val * large_val, large_val * large_val]
    print(f"\nExpected result (Large Numbers):")
    print(f"  {expected}")
    
    if multi_res == expected:
        print("\nMulti-Modulus Verification SUCCESS! ✅")
        print("We successfully reconstructed coefficients > 2^60!")
    else:
        print("\nMulti-Modulus Verification FAILED! ❌")
        
    print("-" * 60)
