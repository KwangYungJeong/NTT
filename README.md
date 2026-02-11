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

This example demonstrates the full cycle of multiplying $A(x)=1+2x+3x^2$ and $B(x)=4+5x$ using $p=17$ and $N=8$.

| Step | Operation | Result / State |
| :--- | :--- | :--- |
| **1. Padding** | Pad to $N=8$ | $A: [1, 2, 3, 0, 0, 0, 0, 0]$, $B: [4, 5, 0, 0, 0, 0, 0, 0]$ |
| **2. F-NTT** | Transform $A, B$ ($w=9$) | $A_{freq}: [6, 12, 12, 8, 2, 7, 7, 5]$, $B_{freq}: [9, 14, 4, 1, 16, 7, 13, 11]$ |
| **3. Multiply** | Point-wise ($A \cdot B$) | $C_{freq}: [3, 15, 14, 8, 15, 15, 6, 4] \pmod{17}$ |
| **4. I-NTT** | Inverse ($w^{-1}=2$) | $C_{coeff}: [4, 13, 5, 15, 0, 0, 0, 0]$ (Note: $22 \equiv 5 \pmod{17}$) |
| **5. Trimming** | Final Result | **$4 + 13x + 5x^2 + 15x^3$** |

> [!IMPORTANT]
> **Numerical Aliasing**: Notice the $x^2$ coefficient is **5**, not 22. This is because $22 \pmod{17} = 5$. For larger results, use a larger prime or **Multi-Modulus NTT**.

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

Solve a system of congruences $x \equiv a_i \pmod{m_i}$ where the moduli are pairwise coprime.

### Run
```bash
python3 src/CRT.py
```

### Methods
-   **Gauss Construction**: Uses the formula $x = \sum a_i M_i y_i \pmod M$, where $M_i = M/m_i$ and $y_i = M_i^{-1} \pmod{m_i}$.
-   **Mixed Radix Conversion (MRC)**: A Garner-style iterative construction $x = v_1 + v_2 m_1 + v_3 m_1 m_2 + \dots$ that is often more efficient for large numbers.

### Worked Example (Side-by-Side)
For Moduli $m = [3, 5, 7, 11]$ and Remainders $a = [2, 3, 2, 4]$:

| Step | Gauss construction ($x = \sum a_i M_i y_i$) | Mixed Radix Conversion (Garner's) |
| :--- | :--- | :--- |
| **1** | $M_1 = 385, y_1 = 1 \Rightarrow w_1 = 2 \cdot 385 \cdot 1 = \mathbf{770}$ | $v_1 = a_1 = \mathbf{2}$ (x_curr: 2) |
| **2** | $M_2 = 231, y_2 = 1 \Rightarrow w_2 = 3 \cdot 231 \cdot 1 = \mathbf{693}$ | $v_2 = (a_2 - x_1) \cdot M_1^{-1} \equiv \mathbf{2} \pmod 5$ (x_curr: 8) |
| **3** | $M_3 = 165, y_3 = 3 \Rightarrow w_3 = 2 \cdot 165 \cdot 3 = \mathbf{990}$ | $v_3 = (a_3 - x_2) \cdot (M_1 m_2)^{-1} \equiv \mathbf{1} \pmod 7$ (x_curr: 23) |
| **4** | $M_4 = 105, y_4 = 2 \Rightarrow w_4 = 4 \cdot 105 \cdot 2 = \mathbf{840}$ | $v_4 = (a_4 - x_3) \cdot (M_1 m_2 m_3)^{-1} \equiv \mathbf{6} \pmod{11}$ (x_curr: 653) |

**Final Result**: Both methods yield $x = \mathbf{653} \pmod{1155}$.

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

---

## References

### NTT & Polynomials
- [Number-theoretic transform (Wiki)](https://en.wikipedia.org/wiki/Number-theoretic_transform)
- [CP-Algorithms: NTT](https://cp-algorithms.com/algebra/fft.html)
- [Project Nayuki: Integer DFT](https://www.nayuki.io/page/number-theoretic-transform-integer-dft)
- [Choosing NTT Primes (Hatsya)](https://hatsya.com/blog/an-efficient-prime-for-number-theoretic-transforms)

### CRT (Chinese Remainder Theorem)
- [CRT Overview (Wiki)](https://en.wikipedia.org/wiki/Chinese_remainder_theorem)
- [Constructive CRT (Gauss-style)](https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Existence_(constructive_proof))
- [Garner’s Algorithm (Mixed Radix)](https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Garner's_algorithm)
- [CP-Algorithms: CRT](https://cp-algorithms.com/algebra/chinese-remainder-theorem.html)

### Mathematical Background
- [Proth Prime (NTT-friendly)](https://en.wikipedia.org/wiki/Proth_prime)
- [Solinas Prime (Goldilocks-type)](https://en.wikipedia.org/wiki/Solinas_prime)
- [2-adic Valuation](https://en.wikipedia.org/wiki/P-adic_order)
- [Primitive Root Modulo n](https://en.wikipedia.org/wiki/Primitive_root_modulo_n)
- [Extended Euclidean Algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm)
- [Hamming Weight (Efficiency)](https://en.wikipedia.org/wiki/Hamming_weight)
