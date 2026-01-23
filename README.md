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

## Worked example (side‑by‑side)

For the default system

- Moduli: \(m = [3, 5, 7, 11]\)
- Remainders: \(a = [2, 3, 2, 4]\)

we can view the Gauss and MRC steps side‑by‑side:

| Step | Gauss construction (using \(x = \sum a_i M_i y_i \pmod M\)) | Mixed Radix Conversion (using \(x = v_1 + v_2 m_1 + v_3 m_1 m_2 + v_4 m_1 m_2 m_3\)) |
|------|-------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| 1 | \(M_1 = 5 \cdot 7 \cdot 11 = 385\)  \newline \(385 \cdot y_1 \equiv 1 \pmod 3 \Rightarrow y_1 = 1\)  \newline \(w_1 = 2 \cdot 385 \cdot 1 = \mathbf{770}\) | \(v_1 = a_1 = \mathbf{2}\)  \newline (current value: \(2\)) |
| 2 | \(M_2 = 3 \cdot 7 \cdot 11 = 231\)  \newline \(231 \cdot y_2 \equiv 1 \pmod 5 \Rightarrow y_2 = 1\)  \newline \(w_2 = 3 \cdot 231 \cdot 1 = \mathbf{693}\) | \(2 + 3 v_2 \equiv 3 \pmod 5\)  \newline \(3 v_2 \equiv 1 \pmod 5 \Rightarrow v_2 = \mathbf{2}\)  \newline (current value: \(2 + 2 \cdot 3 = 8\)) |
| 3 | \(M_3 = 3 \cdot 5 \cdot 11 = 165\)  \newline \(165 \cdot y_3 \equiv 1 \pmod 7 \Rightarrow y_3 = 3\)  \newline \(w_3 = 2 \cdot 165 \cdot 3 = \mathbf{990}\) | \(8 + (3 \cdot 5) v_3 \equiv 2 \pmod 7\)  \newline \(1 + 1 v_3 \equiv 2 \pmod 7 \Rightarrow v_3 = \mathbf{1}\)  \newline (current value: \(8 + 1 \cdot 15 = 23\)) |
| 4 | \(M_4 = 3 \cdot 5 \cdot 7 = 105\)  \newline \(105 \cdot y_4 \equiv 1 \pmod{11} \Rightarrow y_4 = 2\)  \newline \(w_4 = 4 \cdot 105 \cdot 2 = \mathbf{840}\) | \(23 + (3 \cdot 5 \cdot 7) v_4 \equiv 4 \pmod{11}\)  \newline \(1 + 6 v_4 \equiv 4 \pmod{11} \Rightarrow v_4 = \mathbf{6}\)  \newline (current value: \(23 + 6 \cdot 105 = 653\)) |

Final aggregation for each method:

- Gauss: \(x = (770 + 693 + 990 + 840) \pmod{1155} = 3293 \pmod{1155} = \mathbf{653}\)
- MRC: \(x = 2 + 6 + 15 + 630 = \mathbf{653}\)

## References (good starting points)

- Chinese Remainder Theorem (overview): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem`
- Constructive CRT (Gauss-style reconstruction formula): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Existence_(constructive_proof)`
- CRT in competitive programming (constructive method + implementation notes): `https://cp-algorithms.com/algebra/chinese-remainder-theorem.html`
- Garner’s algorithm (mixed radix CRT): `https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Garner's_algorithm`
- Modular inverse / Extended Euclidean algorithm: `https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm`


