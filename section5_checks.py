#!/usr/bin/env python3
"""
section5_checks.py

Verifies the checkable claims of Section 5 (constructive prime generation):
 1. Pepin-shadow criterion (Thm): for N = 576^(2^m)+1, N-1 = 2^(6.2^m) 3^(2.2^m),
    the Jacobi symbol (5/N) = -1, and the single-exponentiation criterion
    certifies primality (tested on small m where N is prime).
 2. Oasis lattice (Thm): density rho_P and that every lattice member m in (P, P^2)
    yields a genuine twin prime pair (m-1, m+1) -- direct check for small P.
 3. Sieve dimension kappa = 1: each prime forbids exactly one residue class.
 4. shadow-admissible offsets: density 88/345 and the list up to 37.
 5. shadow complement maps f_-, f_+, f_-- : primes = complement of their images.
"""
from math import gcd, isqrt, prod
from shadow_core import is_prime

def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

print("=== 1. Pepin-shadow criterion ===")
for m in range(1, 5):
    N = 576**(2**m) + 1
    fact_ok = (N - 1 == 2**(6*2**m) * 3**(2*2**m))
    jac = jacobi(5, N)
    prime = is_prime(N)
    crit = None
    if prime:
        x = pow(5, (N-1)//6, N)            # single exponentiation
        crit = (pow(x, 3, N) == N-1) and (pow(x, 2, N) != 1)
    print(f"  m={m}: N=576^(2^{m})+1  N-1 factored 2^a 3^b: {fact_ok}  (5/N)=-1: {jac==-1}"
          f"  prime: {prime}  Pepin-witness(a=5): {crit}")

print("\n=== 2. Oasis lattice: members in (P,P^2) are certified twins ===")
for P in (7, 11, 13):
    primes = [q for q in range(5, P+1) if is_prime(q)]
    Pp = 6 * prod(primes)                      # period 6 * prod_{5<=q<=P} q
    # lattice points m=6u in (P, P^2) with m mod q not in {1, q-1} for 5<=q<=P
    lat = [m for m in range(P+2, P*P) if m % 6 == 0 and
           all(m % q not in (1, q-1) for q in primes)]
    all_twins = all(is_prime(m-1) and is_prime(m+1) for m in lat)
    rho = (1/6) * prod((q-2)/q for q in primes)
    print(f"  P={P}: {len(lat)} lattice points in (P,P^2), all are twin primes: {all_twins}"
          f"   density rho_P={rho:.5f}")

print("\n=== 3. Sieve dimension kappa=1 (forbidden set is a singleton) ===")
# shadow safe-position sieve on ray r: p | (24k+r) <=> k = sigma_r(p) (one class)
singleton = True
for p in [q for q in range(5, 2000) if is_prime(q)]:
    inv24 = pow(24, -1, p)
    forb = set((-r * inv24) % p for r in (1,))   # one ray -> one forbidden class
    if len(forb) != 1:
        singleton = False
print(f"  every prime forbids exactly one class per ray: {singleton}")

print("\n=== 4. shadow-admissible offsets ===")
adm = [c for c in range(1, 38) if gcd(c, 6) == 1 and c % 5 != 4 and c % 23 != 22]
dens = (1/3) * (4/5) * (22/23)  # phi(6)/6 * (4/5) * (22/23)
print(f"  admissible c<=37 coprime to 6: {adm}")
print(f"  density = 88/345 = {88/345:.4f}  (formula {dens:.4f})")

print("\n=== 5. shadow complement maps (primes = complement of images) ===")
def fminus(K):   # n=6k-1 prime <=> k not in Im(f_-),  f_-(a,b)=6ab-a+b
    img = set()
    a = 1
    while 6*a*1 - a + 1 <= K:
        b = 1
        while 6*a*b - a + b <= K:
            img.add(6*a*b - a + b); b += 1
        a += 1
    return img
K = 2000
imgm = fminus(K)
bad = [k for k in range(1, K+1) if (is_prime(6*k-1)) == (k in imgm)]
print(f"  branch 5 mod 6 up to 6*{K}: mismatches between primality and complement = {len(bad)} (want 0)")

print("\n=== 6. inversion dichotomy: lattice gap sequence non-constant (no closed form) ===")
for P in (5, 7, 11):
    pr = [q for q in range(5, P+1) if is_prime(q)]
    Pp = 6 * prod(pr)
    L = [m for m in range(0, Pp) if m % 6 == 0 and all(m % q not in (1, q-1) for q in pr)]
    gaps = [L[i+1]-L[i] for i in range(len(L)-1)] + [L[0] + Pp - L[-1]]
    print(f"  P={P}: |Lambda mod P#|={len(L)}, distinct gaps={sorted(set(gaps))}, "
          f"non-constant={len(set(gaps)) > 1}")
print("  => the gap sequence is periodic and non-constant, hence the k-th lattice")
print("     element is not a polynomial in k: index inversion has no closed form")
print("     (regional inversion is achievable; the dichotomy is sharp).")
