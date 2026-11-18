"""Numerical sanity checks for the setwise/TV note. Exact rational arithmetic where possible.
Not part of any proof; each check targets one step that could fail."""
from fractions import Fraction as Fr
from itertools import product
import math, random

ok = True
def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)

# 1. Bernoulli example. Y fair two-sided bits conditioned on Y_0 = 1; (T^n Y)_j = Y_{j+n}.
# Law of the window (T^n Y)_{-w..w} versus fair law on {0,1}^{2w+1}.
def window_law(n, w):
    law = {}
    for bits in product((0, 1), repeat=2 * w + 1):
        # coordinate j of T^n Y is Y_{j+n}; forced coordinate is j = -n
        p = Fr(1)
        for idx, j in enumerate(range(-w, w + 1)):
            if j + n == 0:
                p *= 1 if bits[idx] == 1 else 0
            else:
                p *= Fr(1, 2)
        law[bits] = p
    return law

for n in range(0, 6):
    for w in range(0, 5):
        law = window_law(n, w)
        tv = sum(abs(p - Fr(1, 2 ** (2 * w + 1))) for p in law.values()) / 2
        expect = Fr(1, 2) if n <= w else Fr(0)
        check(f"Bernoulli window TV n={n} w={w} is {expect}", tv == expect)

# 2. Criterion (c) on a finite Markov example: Q stationary 3-state Markov chain, h = h(x_0, x_1).
# Then E[h|F_n] = E[h | x_n] for n >= 1 (Markov), so TV_n = (1/2) E|E[h|x_n] - 1|.
# Brute force: sup over events in sigma(x_n..x_{n+L}) of P(theta_n X in A) - Q(A), which is the true TV here
# because E[h|F_n] depends on x_n only.
P = [[Fr(1, 2), Fr(1, 3), Fr(1, 6)], [Fr(1, 4), Fr(1, 4), Fr(1, 2)], [Fr(1, 3), Fr(1, 3), Fr(1, 3)]]
# stationary law by solving pi P = pi
import itertools
def stationary(P):
    # power iteration with fractions converges slowly; solve linear system instead
    a, b, c = P
    # pi = pi P, sum 1 ; solve 3x3 via Cramer
    M = [[P[0][0] - 1, P[1][0], P[2][0]], [P[0][1], P[1][1] - 1, P[2][1]], [1, 1, 1]]
    rhs = [0, 0, 1]
    def det(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    D = det(M)
    pi = []
    for k in range(3):
        Mk = [row[:] for row in M]
        for i in range(3):
            Mk[i][k] = rhs[i]
        pi.append(det(Mk) / D)
    return pi
pi = stationary(P)
check("stationary law sums to 1 and is invariant", sum(pi) == 1 and all(sum(pi[i] * P[i][j] for i in range(3)) == pi[j] for j in range(3)))
h = {(0, 0): Fr(2), (0, 1): Fr(0), (0, 2): Fr(1), (1, 0): Fr(1, 2), (1, 1): Fr(3), (1, 2): Fr(0), (2, 0): Fr(1), (2, 1): Fr(1), (2, 2): Fr(1)}
Eh = sum(pi[i] * P[i][j] * h[(i, j)] for i in range(3) for j in range(3))
h = {k: v / Eh for k, v in h.items()}
def matpow(P, n):
    R = [[Fr(int(i == j)) for j in range(3)] for i in range(3)]
    for _ in range(n):
        R = [[sum(R[i][k] * P[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return R
prev = None
for n in range(1, 6):
    # law of x_n under P=hQ: sum_{i,j} pi_i P_ij h_ij P^{n-1}_{j,s}
    Pn1 = matpow(P, n - 1)
    lawn = [sum(pi[i] * P[i][j] * h[(i, j)] * Pn1[j][s] for i in range(3) for j in range(3)) for s in range(3)]
    # g(s) = E_Q[h | x_n = s] = lawn[s]/pi[s]
    tv_formula = sum(abs(lawn[s] - pi[s]) for s in range(3)) / 2
    # brute force over events in sigma(x_n, x_{n+1}): both laws Markov from x_n with same kernel
    L = 2
    paths = list(product(range(3), repeat=L))
    def prob(start_law, path):
        p = start_law[path[0]]
        for a, b in zip(path, path[1:]):
            p *= P[a][b]
        return p
    diffs = [prob(lawn, q) - prob(pi, q) for q in paths]
    tv_brute = sum(d for d in diffs if d > 0)
    check(f"Markov criterion n={n}: formula TV == brute-force sup", tv_formula == tv_brute)
    if prev is not None:
        check(f"Markov criterion n={n}: TV nonincreasing", tv_formula <= prev)
    prev = tv_formula

# 3. Affinity of 1 + eps sin(2 pi m x) with Lebesgue on [0,1]: independent of m and < 1.
def affinity(eps, m, N=200000):
    s = 0.0
    for i in range(N):
        x = (i + 0.5) / N
        s += math.sqrt(1 + eps * math.sin(2 * math.pi * m * x))
    return s / N
for eps in (0.3, 1.0):
    vals = [affinity(eps, m) for m in (1, 2, 7, 31)]
    check(f"affinity eps={eps}: values {['%.9f' % v for v in vals]} equal and < 1",
          max(vals) - min(vals) < 1e-7 and max(vals) < 1)

# 4. TV <= sqrt(1 - affinity^2) on random finite distributions.
random.seed(1)
good = True
for _ in range(2000):
    k = random.randint(2, 6)
    p = [random.random() for _ in range(k)]; q = [random.random() for _ in range(k)]
    if random.random() < 0.3:
        p[0] = 0.0
    sp, sq = sum(p), sum(q)
    p = [x / sp for x in p]; q = [x / sq for x in q]
    tv = sum(abs(a - b) for a, b in zip(p, q)) / 2
    rho = sum(math.sqrt(a * b) for a, b in zip(p, q))
    good &= tv <= math.sqrt(max(0.0, 1 - rho * rho)) + 1e-12
check("TV <= sqrt(1-rho^2) on 2000 random pairs", good)

# 5. Period-2 chain: X starts at 0 and alternates; X' stationary uniform start.
# P(theta_n X in {x_0 = 0}) alternates 1, 0, while Q(x_0=0) = 1/2; Cesaro TV -> 0.
for N in (1, 2, 3, 10, 11):
    # laws of theta_k X are point masses on the two alternating paths a (start 0) and b (start 1)
    avgX = {'a': Fr(sum(1 for k in range(N) if k % 2 == 0), N), 'b': Fr(sum(1 for k in range(N) if k % 2 == 1), N)}
    tv = (abs(avgX['a'] - Fr(1, 2)) + abs(avgX['b'] - Fr(1, 2))) / 2
    expect = Fr(0) if N % 2 == 0 else Fr(1, 2 * N)
    check(f"period-2 Cesaro TV N={N} equals {expect}", tv == expect)

print('ALL PASS' if ok else 'SOME FAIL')
