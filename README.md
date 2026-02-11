# NTT & Foundation Tools

A high-performance implementation of the **Number Theoretic Transform (NTT)** and essential modular arithmetic tools, including **CRT** (Chinese Remainder Theorem) and **Prime Search** utilities.

## Requirements
- Python 3.x
- `sympy` library (for primality testing and primitive root discovery)
  ```bash
  pip install sympy
  ```

---

## 1. Quick Start: NTT Demo

The core of this project is a robust $O(N \log N)$ polynomial multiplication engine.

### Run
```bash
# Run with default prime (469762049)
python3 src/NTT.py

# Run with a custom prime (e.g., 2013265921)
# The primitive root (g) is calculated automatically.
python3 src/NTT.py --prime 2013265921
```

### Dynamic Configuration
The system **automatically discovers** the required mathematical parameters:
- **Automatic Root Detection**: Primitive root $g$ is calculated via `sympy` in real-time.
- **Safety Checks**: Automatically verifies if the prime's **2-adic valuation** supports the target polynomial size $N$.

---

## 2. NTT Step-by-Step Breakdown

NTT enables fast multiplication by switching between the **Coefficient Domain** and the **Frequency Domain**.

### Steps 1-5
1.  **Size & Padding**: Find smallest power of 2 ($N$) and pad with zeros to avoid Aliasing.
2.  **Forward NTT**: Transform using: $X_k = \sum_{j=0}^{N-1} a_j \cdot \omega^{kj} \pmod p$
3.  **Point-wise Multiply**: Calculate $C_{freq}[i] = (A_{freq}[i] \cdot B_{freq}[i]) \pmod p$ in $O(N)$.
4.  **Inverse NTT**: Restore coefficients: $a_j = N^{-1} \sum_{k=0}^{N-1} X_k \cdot \omega^{-kj} \pmod p$
5.  **Trimming**: Remove trailing zeros to get the final result.

### Worked Example ($N=8, p=17$)
| Stage | Process | Data State (Bit-reversed) |
| :--- | :--- | :--- |
| **Initial** | Input coefficients | `[1, 2, 3, 4, 0, 0, 0, 0]` |
| **Bit-Reversal** | Reorder indices | `[1, 0, 3, 0, 2, 0, 4, 0]` |
| **F-NTT** | Butterfly with $w=9$ | **`[10, 16, 6, 11, 15, 13, 7, 15]`** |

---

## 3. Large Integer Multiplication (Multi-Modulus NTT)

When coefficients exceed a single modulus (e.g., $> 2^{31}$), we perform NTT over multiple primes and reconstruct the final result using **CRT**.

### Run
```bash
python3 src/multi_mod_ntt.py
```
-   **Capability**: Successfully reconstructs coefficients $> 2^{60}$ using 3-moduli sets.

---

## 4. Prime Search Tools

Find optimal primes for your specific NTT or cryptographic needs.

### NTT Primes
Search for primes of the form $k \cdot 2^n + 1$:
```bash
python3 src/prime_search.py search_ntt_prime --n_power 20 --count 5
```

### Goldilocks/Solinas Primes
Find primes like $2^n - c$ for fast modular reduction:
```bash
python3 src/prime_search.py search_goldilock_prime --n_start 20 --n_end 31
```

---

## 5. CRT (Chinese Remainder Theorem) Demo

Independent utility to solve systems of congruences.

### Run
```bash
python3 src/CRT.py
```
-   **Methods**: Includes both **Gauss construction** and **Mixed Radix Conversion (MRC)**.

---

## 6. Technical Constraints & Design

### NTT Constraints
| Constraint | Decided By | Consequence if Exceeded |
| :--- | :--- | :--- |
| **Horizontal (N)** | 2-adic valuation | Calculation Error (Root doesn't exist) |
| **Vertical ($C_i$)** | Modulus Magnitude | Numerical Aliasing (Wrap-around) |

### Choosing a Good NTT Prime
1.  **Form**: $p = k \cdot 2^n + 1$ (Proth prime).
2.  **Magnitude**: Must be larger than the max coefficient (unless using Multi-Mod).
3.  **Efficiency**: Fermat numbers or Low-popcount primes allow for faster reduction.

---

## Project Structure
- `src/NTT.py`: Core NTT engine with dynamic config.
- `src/multi_mod_ntt.py`: Orchestrates multi-prime multiplication.
- `src/CRT.py`: CRT solvers (Gauss, MRC).
- `src/prime_search.py`: Math utilities and prime discovery.

---

## References
- [Number-theoretic transform (Wiki)](https://en.wikipedia.org/wiki/Number-theoretic_transform)
- [CP-Algorithms: NTT](https://cp-algorithms.com/algebra/fft.html)
- [Project Nayuki: Integer DFT](https://www.nayuki.io/page/number-theoretic-transform-integer-dft)
