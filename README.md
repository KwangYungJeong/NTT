# NTT & Foundation Tools

## Requirements

- Python 3.x (no third-party dependencies)
- `sympy` (for primality testing and modular inverse)


## CRT (Chinese Remainder Theorem) Demo

This repo contains a single Python script, `CRT.py`, that demonstrates solving a system of congruences using:

- **Gauss construction** (classic CRT reconstruction via $\sum a_i M_i y_i \bmod M$)
- **Mixed Radix / MRC** (a Garner-style step-by-step construction of the solution)

The script prints intermediate steps so you can follow the arithmetic.

### Run

```bash
python3 src/CRT.py
```

By default, the script uses:

- Moduli: `[3, 5, 7, 11]`
- Remainders: `[2, 3, 2, 4]`

and prints the final solution $x$ (mod $M$) computed by both methods.

### Notes

- CRT requires the moduli to be **pairwise coprime** for a unique solution modulo $M = \prod m_i$.
- The MRC section maintains a running partial solution `x_curr` and updates it at each step; this makes the final “reconstruction” step unnecessary.

### Worked example (side‑by‑side)

For the default system

- Moduli: $m = [3, 5, 7, 11]$
- Remainders: $a = [2, 3, 2, 4]$

we can view the Gauss and MRC steps side‑by‑side:

| Step | Gauss construction (using $x = \sum a_i M_i y_i \pmod M$) | Mixed Radix Conversion (using $x = v_1 + v_2 m_1 + v_3 m_1 m_2 + v_4 m_1 m_2 m_3$) |
|------|-----------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| 1 | $M_1 = 5 \cdot 7 \cdot 11 = 385$<br>$385 \cdot y_1 \equiv 1 \pmod 3 \Rightarrow y_1 = 1$<br>$w_1 = 2 \cdot 385 \cdot 1 = \mathbf{770}$ | $v_1 = a_1 = \mathbf{2}$<br>(current value: $2$) |
| 2 | $M_2 = 3 \cdot 7 \cdot 11 = 231$<br>$231 \cdot y_2 \equiv 1 \pmod 5 \Rightarrow y_2 = 1$<br>$w_2 = 3 \cdot 231 \cdot 1 = \mathbf{693}$ | $2 + 3 v_2 \equiv 3 \pmod 5$<br>$3 v_2 \equiv 1 \pmod 5 \Rightarrow v_2 = \mathbf{2}$<br>(current value: $2 + 2 \cdot 3 = 8$) |
| 3 | $M_3 = 3 \cdot 5 \cdot 11 = 165$<br>$165 \cdot y_3 \equiv 1 \pmod 7 \Rightarrow y_3 = 3$<br>$w_3 = 2 \cdot 165 \cdot 3 = \mathbf{990}$ | $8 + (3 \cdot 5) v_3 \equiv 2 \pmod 7$<br>$1 + 1 v_3 \equiv 2 \pmod 7 \Rightarrow v_3 = \mathbf{1}$<br>(current value: $8 + 1 \cdot 15 = 23$) |
| 4 | $M_4 = 3 \cdot 5 \cdot 7 = 105$<br>$105 \cdot y_4 \equiv 1 \pmod{11} \Rightarrow y_4 = 2$<br>$w_4 = 4 \cdot 105 \cdot 2 = \mathbf{840}$ | $23 + (3 \cdot 5 \cdot 7) v_4 \equiv 4 \pmod{11}$<br>$1 + 6 v_4 \equiv 4 \pmod{11} \Rightarrow v_4 = \mathbf{6}$<br>(current value: $23 + 6 \cdot 105 = 653$) |

Final aggregation for each method:

- Gauss: $x = (770 + 693 + 990 + 840) \pmod{1155} = 3293 \pmod{1155} = \mathbf{653}$
- MRC: $x = 2 + 6 + 15 + 630 = \mathbf{653}$

### References (good starting points)

- Chinese Remainder Theorem (overview): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem`
- Constructive CRT (Gauss-style reconstruction formula): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Existence_(constructive_proof)`
- CRT in competitive programming (constructive method + implementation notes): `https://cp-algorithms.com/algebra/chinese-remainder-theorem.html`
- Garner’s algorithm (mixed radix CRT): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Garner's_algorithm`
- Modular inverse / Extended Euclidean algorithm: `https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm`



## NTT (Number Theoretic Transform) Demo

This script demonstrates the core NTT algorithm ($O(N \log N)$) to multiply two polynomials. It includes detailed console output for each step of the calculation.

### Run

```bash
# Run with default prime (469762049)
python3 src/NTT.py

# Run with a custom prime (e.g., 17)
# Primitive root (g) will be calculated automatically
python3 src/NTT.py --prime 17
```

### Default Parameters
By default, the script uses a production-grade **NTT-friendly prime**:

- **Modulus ($p$)**: `469762049` ($7 \cdot 2^{26} + 1$)
- **Primitive Root ($g$)**: `3`
- **Capacity ($N$)**: Up to $2^{26}$ (~67 million coefficients).

---

## Prime Search Tools

This project includes a Python script `prime_search.py` to find suitable primes for NTT and other modular arithmetic applications.

### 1. Search for NTT Primes
Finds primes of the form $p = c \cdot 2^k + 1$ that support large NTT sizes.

```bash
# Find 5 primes with at least 2^20 capacity
python3 src/prime_search.py search_ntt_prime --n_power 20 --count 5
```

**Key Output:**
- **2-adic valuation**: Max power of 2 that divides $p-1$. Determines max FFT size.
- **Max Root of Unity ($w$)**: Primitive root of unity for the max FFT size.

### 2. Search for Goldilocks/Solinas Primes
Finds primes of the form $2^n - c$ or $2^n - 2^m \pm 1$ optimized for fast modular reduction (but usually bad for NTT).

```bash
# Search for Solinas primes near 2^20 to 2^31
python3 src/prime_search.py search_goldilock_prime --n_start 20 --n_end 31
```

**Key Output:**
- **Golden! (phi)**: Indicates $2m \approx n$, allowing for Karatsuba-like optimization.
- **Word 32-bit**: Indicates $n$ is a multiple of 32, suitable for word-aligned implementation.

---

## NTT Step-by-Step Breakdown

The Number Theoretic Transform (NTT) enables fast polynomial multiplication by switching between the **Coefficient Domain** and the **Frequency (Point-Value) Domain**.

### Step 1: Size Determination & Padding
1.  Calculate target length: $L = len(A) + len(B) - 1$.
2.  Find $N$ as the smallest power of 2 such that $N \ge L$.
3.  Pad both polynomials with zeros up to length $N$.
    -   This prevents **Cyclic Convolution (Aliasing)**.

### Step 2: Forward NTT (Transform)
Moves polynomials to the Frequency Domain.
1.  **Bit-Reversal Permutation**: Reorders indices for iterative processing.
2.  **Butterfly Operations**: Combines values using the **summation formula**:
    $$X_k = \sum_{j=0}^{N-1} a_j \cdot \omega^{kj} \pmod p$$
    (where $X_k$ is the evaluation at $\omega^k$)

### Step 3: Point-wise Multiplication
In the Frequency Domain, multiplication of two polynomials is simply the element-by-element product of their evaluations:
-   $C_{freq}[i] = (A_{freq}[i] \cdot B_{freq}[i]) \pmod p$
-   Complexity: **$O(N)$** (Much faster than $O(N^2)$ convolution!)

### Step 4: Inverse NTT (Reverse Transform)
Moves the product back to the Coefficient Domain utilizing the **inverse formula**:
$$a_j = N^{-1} \sum_{k=0}^{N-1} X_k \cdot \omega^{-kj} \pmod p$$

1.  Uses the same butterfly structure but with the **inverse root of unity** $\omega^{-1}$.
2.  Requires a final **Scaling**: multiply all elements by $N^{-1} \pmod p$.

### Step 5: Trimming
Remove the trailing zeros beyond the target length $L$ to get the final coefficients.

---

### Worked Example: $N=8, p=17$
- **Input**: $A(x) = 1 + 2x + 3x^2 + 4x^3$
- **Parameters**: Modulo $p=17$, 8th root of unity $w=9$, $w^{-1}=2$.

| Stage | Process | Data State (Indices bit-reversed) |
| :--- | :--- | :--- |
| **Initial** | Input coefficients | `[1, 2, 3, 4, 0, 0, 0, 0]` |
| **Bit-Reversal** | Reorder indices | `[1, 0, 3, 0, 2, 0, 4, 0]` |
| **F-NTT (Len 8)** | Butterfly with $w=9$ | **`[10, 16, 6, 11, 15, 13, 7, 15]`** |

**Verification (Frequency Domain)**: 
The NTT result $X_k$ corresponds to the evaluation $A(w^k) \pmod p$.
- $A(1) = 1+2+3+4 = \mathbf{10}$
- $A(9) = 1(9)^0 + 2(9)^1 + 3(9)^2 + 4(9)^3 \equiv \mathbf{16}$
- $A(13) = 1(13)^0 + 2(13)^1 + 3(13)^2 + 4(13)^3 \equiv \mathbf{6}$

---

## NTT Constraints & Limits

When designing an NTT-based system, you must consider two primary mathematical constraints determined by your choice of modulus $m$.

### 1. Horizontal Limit: Number of Coefficients ($N$)
- **Constraint**: The maximum transform size $N$ is limited by the **2-adic valuation** of $m-1$.
- **Formula**: if $m-1 = k \cdot 2^n$, then $N_{max} = 2^n$.
- **In this Repo**: $m = 469762049 = 7 \cdot 2^{26} + 1$, so $N \le 2^{26}$ (~67 million).
- **Consequence**: Attempting a transform larger than this will fail as the required $N$-th roots of unity do not exist.

### 2. Vertical Limit: Magnitude of Coefficients ($C_i$)
- **Constraint**: Each resulting coefficient $C_i$ must be smaller than the modulus $m$ to avoid **modular wrap-around (aliasing)**.
- **Formula**: $\max(C_i) < m$. 
- **In this Repo**: Results $> 469,762,049$ will be "cut" by the modulo (e.g., $5 \cdot 10^8 \pmod m = 30,237,951$).
- **Solution**: If your expected coefficients are larger than $m$, you must use **Multi-Modulus NTT** (combine results from $m_1, m_2, \dots$ via CRT).

| Constraint | Decided By | Default Limit | If Exceeded? |
| :--- | :--- | :--- | :--- |
| **Data Size (N)** | $m-1$ structure | $2^{26}$ elements | Calculation Error |
| **Value ($C_i$)** | $|m|$ magnitude | $469,762,049$ | Numerical Aliasing |

---

---

## Choosing a Good NTT Prime

A "good" prime for NTT must satisfy specific mathematical and computational criteria:

1.  **Form: $p = k \cdot N + 1$**: This guarantees the existence of an $N$-th primitive root of unity. For FFT-like efficiency, $N$ should be a power of 2 ($2^n$).
2.  **Magnitude**: $p$ must be larger than the maximum possible coefficient in the result polynomial to avoid aliasing.
3.  **Efficiency**: Primes like **Proth primes** or **Fermat numbers** ($2^k + 1$) often allow for faster modular reduction.

### Key References (Deep Dive)
- **Project Nayuki**: [Number-theoretic transform (integer DFT)](https://www.nayuki.io/page/number-theoretic-transform-integer-dft) - A great overview of selection criteria.
- **Kyber/Dilithium NTT**: Primes like `12289` (BLISS) and `3329` (Kyber) are standards in PQC.
- **Choosing NTT Primes**: [How to find good primes](https://hatsya.com/blog/an-efficient-prime-for-number-theoretic-transforms) - Technical blog on large NTT primes.
- **Arithmetic in Finite Fields**: Standard textbooks (e.g., *Menezes et al.*) cover primitive roots and modular arithmetic.

---

### References (Prime & NTT Basics)

- **Number Theoretic Transform (NTT)**: `https://en.wikipedia.org/wiki/Number-theoretic_transform`
- **Proth Prime** (NTT-friendly): `https://en.wikipedia.org/wiki/Proth_prime`
- **Solinas Prime** (Goldilocks-type): `https://en.wikipedia.org/wiki/Solinas_prime`
- **2-adic valuation**: `https://en.wikipedia.org/wiki/P-adic_order`
- **Primitive Root**: `https://en.wikipedia.org/wiki/Primitive_root_modulo_n`
