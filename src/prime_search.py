
import math
import argparse
import sys

try:
    from sympy import isprime
    from sympy.ntheory import primitive_root
except ImportError:
    print("Error: 'sympy' library is not installed.")
    print("Please install it using: pip install sympy")
    sys.exit(1)

def search_ntt_prime(n_power=20, count=5, lower_g=False):
    """
    Finds primes suitable for Number Theoretic Transform (NTT).
    Form: p = k * 2^n_power + 1
    """
    chk_val = 1 << n_power
    print(f"Searching for NTT primes of the form p = k * 2^{n_power} + 1 (N >= {chk_val})")
    print(f"Start searching... (Max limit: 2^31 - 1)")
    
    candidates = []
    
    k = 1
    while True:
        p = k * chk_val + 1
        
        # Stop condition: p must fit in 32-bit signed integer
        if p >= 2**31:
            break
            
        if isprime(p):
            is_optimal = True
            reasons = []

            # Condition: Check exact power of 2 divisibility
            p_minus_1 = p - 1
            exact_k = 0
            while p_minus_1 % 2 == 0:
                p_minus_1 //= 2
                exact_k += 1
            
            # Find a primitive root modulo p
            g = primitive_root(p)
            
            # Condition: Small primitive root (optional but good)
            if lower_g and g > 10:
                is_optimal = False
                reasons.append(f"g={g} > 10")

            # Condition: Low Hamming Weight (Popcount)
            popcount = bin(p).count('1')
            
            if is_optimal:
                score = (exact_k * 1000) - (popcount * 10) - g
                candidates.append({
                    "p": p, 
                    "exact_k": exact_k, 
                    "popcount": popcount, 
                    "g": g,
                    "score": score
                })
        k += 1

    # Post-processing
    print("\n" + "="*60)
    print(f"Scanning Complete! Found {len(candidates)} candidates in 32-bit range.")
    print(f"Selecting Top {count} Primes...")
    print("Sorting Criteria: Max 2^k > Low Popcount > Small g")
    print("="*60)

    candidates.sort(key=lambda x: x["score"], reverse=True)
    selected = candidates[:count]
    
    for idx, cand in enumerate(selected):
        p = cand["p"]
        max_N = 1 << cand['exact_k']
        w_N = pow(cand['g'], (p - 1) // max_N, p)
        
        print(f"Rank {idx+1}: {p}")
        print(f"    - 2-adic valuation: 2^{cand['exact_k']} (Max N: {max_N})")
        print(f"    - Primitive Root (g): {cand['g']}")
        print(f"    - Max Root of Unity (w_{max_N}): {w_N}")
        print(f"    - Popcount: {cand['popcount']}")
        print(f"    - Hex: {hex(p)}")
        print("-" * 40)

def search_goldilock_prime(n_start=20, n_end=31):
    """
    Finds "Goldilocks-like" primes (Solinas primes) of the form:
    1. 2^n - c (where c is small)
    2. 2^n - 2^m - 1 (Generalized Mersenne)
    """
    print(f"Searching for Goldilocks-like primes near 2^{n_start} to 2^{n_end}...")
    print(f"{'Prime (p)':<15} | {'Form':<20} | {'Type':<15} | {'Note'}")
    print("-" * 80)
    
    candidates = []

    for n in range(n_start, n_end + 1):
        # Type 1: Pseudo-Mersenne (2^n - c)
        for c in range(1, 100):
            p = (1 << n) - c
            if isprime(p):
                note = ""
                if n % 32 == 0 or n % 64 == 0:
                     note += "Word-Aligned "
                
                candidates.append({
                    "p": p,
                    "form": f"2^{n} - {c}",
                    "type": "Pseudo-Mersenne",
                    "note": note,
                    "score": 0 - c # Prefer smaller c
                })
                break 

        # Type 2: Solinas (2^n - 2^m +/- 1)
        for m in range(1, n // 2 + 1):
            # Check both -1 and +1
            for sign in [-1, 1]:
                p = (1 << n) - (1 << m) + sign
                
                # Check 2^31 limit (exclude if >= 2^31)
                if p >= 2**31:
                    continue
                    
                if isprime(p):
                    score = 0
                    note = []
                    
                    # Check for "Golden" ratio property (m approx n/2)
                    if 2 * m == n:
                        note.append("Golden! (phi)")
                        score += 100
                    elif abs(2 * m - n) <= 1:
                         note.append("Near-Golden")
                         score += 50
                    
                    # Check for Word Alignment
                    if n % 32 == 0:
                        note.append("Word 32-bit")
                        score += 20
                    
                    sign_str = "-" if sign == -1 else "+"
                    
                    # Calculate 2-adic valuation
                    p_minus_1 = p - 1
                    k = 0
                    while p_minus_1 % 2 == 0:
                        p_minus_1 //= 2
                        k += 1

                    candidates.append({
                        "p": p,
                        "form": f"2^{n} - 2^{m} {sign_str} 1",
                        "type": "Solinas",
                        "note": ", ".join(note),
                        "score": score,
                        "k": k # Added k to dict for verification
                    })

    # Sort by Score (Desc) then by Prime size (Asc)
    candidates.sort(key=lambda x: (x.get("score", 0), x["p"]), reverse=True)

    for cand in candidates:
        print(f"{cand['p']:<15} | {cand['form']:<20} | {cand['type']:<15} | {cand['note']}")

    print("-" * 80)
    print("\n[Guide]")
    print("1. 'Golden! (phi)': 2m = n. Best for Karatsuba-like multiplication.")
    print("   -> e.g. 2^448 - 2^224 - 1 (Ed448)")
    print("2. 'Word 32-bit': n is a multiple of 32. Fits implementation perfectly.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find primes suitable for NTT/FFT or Goldilocks forms.")
    subparsers = parser.add_subparsers(dest="command", help="Search type")
    
    # Subcommand for NTT primes
    parser_ntt = subparsers.add_parser("search_ntt_prime", help="Search for NTT-friendly primes")
    parser_ntt.add_argument("--n_power", type=int, default=20, help="Minimum power of 2 for FFT size")
    parser_ntt.add_argument("--count", type=int, default=5, help="Number of primes to display")
    parser_ntt.add_argument("--lower_g", action="store_true", help="Prefer smaller primitive roots")
    
    # Subcommand for Goldilocks primes
    parser_gold = subparsers.add_parser("search_goldilock_prime", help="Search for Goldilocks/Solinas primes")
    parser_gold.add_argument("--n_start", type=int, default=20, help="Start power of 2")
    parser_gold.add_argument("--n_end", type=int, default=31, help="End power of 2")
    
    args = parser.parse_args()
    
    if args.command == "search_ntt_prime":
        search_ntt_prime(args.n_power, args.count, args.lower_g)
    elif args.command == "search_goldilock_prime":
        search_goldilock_prime(args.n_start, args.n_end)
    else:
        parser.print_help()
