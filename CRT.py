import math

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return gcd, x, y

# Modular Inverse using Extended Euclidean Algorithm
# a^-1 mod m
# example: mod_inverse(3, 11) = 4 because (3*4) % 11 = 1
# example: mod_inverse(4, 7) = 2 because (4*2) % 7 = 1
# It returns x such that (a*x) % m = 1
# It is essential for CRT calculations for Gauss and MRC methods
def mod_inverse(a, m):
    gcd, x, y = extended_gcd(a, m)
    if gcd != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return (x % m + m) % m


# 1. Gauss's Construction Method
def gauss_crt(m, a):
    total_sum = 0
    M = math.prod(m) # Product of all moduli
    print("--- Gauss Method Steps ---")
    for i in range(len(m)):
        print(f"Step {i+1}:")
        print(f"\tModulus m{i+1} = {m[i]}, Remainder a{i+1} = {a[i]}")
        Mi = M // m[i] ## M_i for each modulus
        print(f"\tCalculating M_i = M / m{i+1} = {M} / {m[i]} = {Mi}")
        yi = mod_inverse(Mi, m[i])
        print(f"\tCalculating y{i+1} = M_i^(-1) mod m{i+1} = {Mi}^(-1) mod {m[i]} = {yi}, because ({Mi} * {yi}) % {m[i]} = 1")
        wi = a[i] * Mi * yi
        print(f"\t==> w{i+1} = {a[i]} * {Mi} * {yi} = {wi}")
        total_sum += wi

    print(f"Final Step: x ≡ Σ w_i (mod M) = {total_sum} (mod {M}) = {total_sum % M}")
    # print(f"Thus, x ≡ {total_sum % M} (mod {M})")
    return total_sum % M


# 2. Mixed Radix Conversion (MRC Method)
def mrc_crt(m, a):
    n = len(m)
    v = [0] * n
    v[0] = a[0]
    print("\n--- MRC Method Steps ---")
    print(f"Step 1:")
    print(f"\tModulus m1 = {m[0]}, Remainder a1 = {a[0]}")
    print(f"\t==> v1 = a1 = {v[0]}")
    print(f"\tPartial solution x1 = v1 = {v[0]}  (mod {m[0]})")
    x_curr = v[0]          # current partial solution x_i
    Mi = m[0]              # product m1*...*m_i used as multiplier in the next step

    for i in range(1, n):
        # At step (i+1), we already have partial solution x_i in x_curr.
        # We use Mi = m1*...*m_i and solve:
        #   x_curr + v_{i+1} * Mi ≡ a_{i+1} (mod m_{i+1})
        x_prev = x_curr

        inv = mod_inverse(Mi, m[i])
        # delta is the difference between the new target remainder a_{i+1}
        # and the current partial solution x_i; we correct x_i by adding
        # v_{i+1} * Mi so that the new congruence is satisfied.
        delta = a[i] - x_prev
        v[i] = (delta * inv) % m[i]
        x_curr = x_prev + v[i] * Mi

        print(f"Step {i+1}:")
        print(f"\tModulus m{i+1} = {m[i]}, Remainder a{i+1} = {a[i]}")
        print(f"\tCurrent partial x{i} = {x_prev}")
        print(f"\tCalculating M{i} = m1*...*m{i} = {Mi}")
        print(
            f"\tCalculating inv = M{i}^(-1) mod m{i+1} = {Mi}^(-1) mod {m[i]} = {inv}, "
            f"because ({Mi} * {inv}) % {m[i]} = 1"
        )
        print(f"\tDelta = a{i+1} - x{i} = {a[i]} - {x_prev} = {delta}  (how far current x is from the new remainder)")
        print(f"\t==> v{i+1} = Delta * inv (mod m{i+1}) = {delta} * {inv} (mod {m[i]}) = {v[i]}")
        print(f"\tUpdate partial x{i+1} = x{i} + v{i+1}*M{i} = {x_prev} + {v[i]}*{Mi} = {x_curr}")

        # Update Mi for the next step: Mi = m1*...*m_{i+1}
        Mi *= m[i]

    print("\n--- MRC Final Result ---")
    print(f"Final x from MRC (mixed-radix) = {x_curr}")
    return x_curr

if __name__ == "__main__":
    # 입력 데이터
    m = [3, 5, 7, 11]
    a = [2, 3, 2, 4]
    # M = math.prod(m)
    print(f"Moduli: {m}")
    print(f"Remainders: {a}")
    print(f"Total Modulus M: {math.prod(m)}")

    # 실행 및 결과 확인
    print("\n--- CRT Calculation ---")
    print("*" * 30)
    print("\nCalculating using Gauss Method:")
    gauss_res = gauss_crt(m, a)
    print("*" * 30)
    print("\nCalculating using MRC Method:")
    mrc_res = mrc_crt(m, a)

    print("*" * 30)
    print(f"\n[Final Result]")
    print(f"Gauss Method: {gauss_res}")
    print(f"MRC Method: {mrc_res}")
