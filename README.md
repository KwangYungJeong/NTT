# CRT (Chinese Remainder Theorem) Demo

This repo contains a single Python script, `CRT.py`, that demonstrates solving a system of congruences using:

- **Gauss construction** (classic CRT reconstruction via \(\sum a_i M_i y_i \bmod M\))
- **Mixed Radix / MRC** (a Garner-style step-by-step construction of the solution)

The script prints intermediate steps so you can follow the arithmetic.

## Requirements

- Python 3.x (no third-party dependencies)

## Run

```bash
python3 CRT.py
```

By default, the script uses:

- Moduli: `[3, 5, 7, 11]`
- Remainders: `[2, 3, 2, 4]`

and prints the final solution \(x\) (mod \(M\)) computed by both methods.

## Notes

- CRT requires the moduli to be **pairwise coprime** for a unique solution modulo \(M = \prod m_i\).
- The MRC section maintains a running partial solution `x_curr` and updates it at each step; this makes the final “reconstruction” step unnecessary.

## References (good starting points)

- Chinese Remainder Theorem (overview): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem`
- Constructive CRT (Gauss-style reconstruction formula): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Existence_(constructive_proof)`
- CRT in competitive programming (constructive method + implementation notes): `https://cp-algorithms.com/algebra/chinese-remainder-theorem.html`
- Garner’s algorithm (mixed radix CRT): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Garner's_algorithm`
- Modular inverse / Extended Euclidean algorithm: `https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm`


