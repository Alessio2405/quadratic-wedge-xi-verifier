"""Actual whitened perturbation ||S^{-1}(M_xi - M_ref)S^{-T}||_2 on real xi coefficients:
3-point (Michalowski log-quadratic) reference vs 4-point q-Pochhammer reference."""
import sys
import mpmath as mp
from xi_data import load

mp.mp.dps = 110
lam, _ = load()


def qp(a, qq, n):
    if n >= 0:
        p = mp.mpf(1)
        for j in range(n):
            p *= 1 - a * qq**j
        return p
    p = mp.mpf(1)
    for j in range(1, -n + 1):
        p *= 1 - a * qq**(-j)
    return 1 / p


def ldl(M):
    n = M.rows
    L = mp.eye(n); D = [0] * n; A = M.copy()
    for k in range(n):
        D[k] = A[k, k]
        for i in range(k + 1, n):
            L[i, k] = A[i, k] / D[k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i, j] -= L[i, k] * D[k] * L[j, k]
    return L, D


def norm2(M):
    ev = mp.eigsy(M.T * M, eigvals_only=True)
    return mp.sqrt(max(abs(e) for e in ev))


def experiment(k, r):
    d = r - 1
    rho = lam[k + 1] / lam[k]
    b = {s: lam[k + s] / (lam[k] * rho**s) for s in range(-max(d, 1), max(d, 2) + 1)}
    u, v = b[-1], b[2]
    tau = -mp.log(u)
    # 4-point q-Pochhammer reference
    qq = (1 - v) / (v * (1 - u)); al = (1 - u) * (1 - v) / (v - u)
    c4 = {s: (1 - al)**s / qp(al, qq, s) for s in range(-max(d, 1), max(d, 2) + 1)}
    for s in (-1, 0, 1, 2):
        assert abs(c4[s] / b[s] - 1) < mp.mpf(10)**-90
    # 3-point log-quadratic reference (Michalowski): c_s = q^{s(s-1)/2}
    c3 = {s: u**(mp.mpf(s * (s - 1)) / 2) for s in range(-d, d + 1)}
    out = []
    for c in (c3, c4):
        M = mp.matrix(r, r); E = mp.matrix(r, r)
        for i in range(r):
            for j in range(r):
                M[i, j] = c[d - i - j]
                E[i, j] = b[d - i - j] - c[d - i - j]
        L, D = ldl(M)
        assert [int(mp.sign(x)) for x in D] == [(-1)**m for m in range(r)]
        S = L * mp.diag([mp.sqrt(abs(x)) for x in D]); Si = S**-1
        out.append(norm2(Si * E * Si.T))
    return out, al, mp.log(qq)


if __name__ == '__main__':
    print('   k    r    ||E||_3pt        ||E||_4pt     ratio3/4   r^2/k        4pt*k/r^2')
    for k in [3000, 6000, 11990]:
        for r in [2, 3, 4, 6, 8, 12, 16]:
            if k + r > 12016:
                continue
            (e3, e4), al, th = experiment(k, r)
            print('%6d %3d   %.4e     %.4e    %8.1f   %.3e    %.4f' % (k, r, e3, e4, e3 / e4, r * r / k, e4 * k / r**2))
