"""Direct numerical test of the transferred Gate B response bound for the q-Pochhammer reference.

Builds M_ref_{ij} = c_{d-i-j} directly from c_s = (1-alpha)^s/(alpha;qq)_s (no formulas reused),
computes S from a *numerical* LDL of M_ref, and measures
    R(s^n) = || S^{-1} (M_ref o s^n) S^{-T} ||_2
against w^n with w = 4 sqrt(r/theta) (claimed) and against the sharper t_true/theta.
"""
import mpmath as mp
import random

mp.mp.dps = 260


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


def c(s, al, qq):
    return (1 - al)**s / qp(al, qq, s)


def ldl(M):
    n = M.rows
    L = mp.eye(n)
    D = [mp.mpf(0)] * n
    A = M.copy()
    for k in range(n):
        D[k] = A[k, k]
        for i in range(k + 1, n):
            L[i, k] = A[i, k] / D[k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i, j] -= L[i, k] * D[k] * L[j, k]
    return L, D


def norm2(M):
    # spectral norm via eigenvalues of M^T M
    ev = mp.eigsy(M.T * M, eigvals_only=True)
    return mp.sqrt(max(abs(e) for e in ev))


def logm_series(X):
    n = X.rows
    N = X - mp.eye(n)
    out = mp.zeros(n)
    P = mp.eye(n)
    for k in range(1, 400):
        P = P * N
        out += ((-1)**(k + 1) / mp.mpf(k)) * P
        if mp.mnorm(P, 1) < mp.mpf(10)**(-80):
            break
    return out


if __name__ == '__main__':
    random.seed(1)
    worst = 0
    rows = []
    for trial in range(12):
        al = mp.mpf(random.uniform(1.34, 2.66))
        theta = mp.mpf(10)**(-random.uniform(2.6, 4.0))
        qq = mp.e**theta
        r = random.choice([3, 5, 8, 10, 12])
        d = r - 1
        if r * theta > mp.mpf(1) / 64:
            continue
        M = mp.matrix(r, r)
        for i in range(r):
            for j in range(r):
                M[i, j] = c(d - i - j, al, qq)
        L, D = ldl(M)
        # sign pattern
        signs = [int(mp.sign(x)) for x in D]
        assert signs == [(-1)**m for m in range(r)], signs
        S = L * mp.diag([mp.sqrt(abs(x)) for x in D])
        Si = S**-1
        # sharper t: theta*d + 2||G||,  G = log Rhat_1, Rhat_1 = S^{-1} diag(Q^i) S  (A-diag commutes)
        Qv = 1 / qq
        Rhat = Si * mp.diag([Qv**i for i in range(r)]) * S
        # Note: S includes c_d*A implicitly (numerical LDL of full M_ref); this checks that A drops out.
        G = logm_series(Rhat)
        t_true = theta * d + 2 * norm2(G)
        t_claim = 4 * mp.sqrt(r * theta)
        assert t_true <= t_claim
        wc = t_claim / theta
        for n in range(0, 9):
            Mn = mp.matrix(r, r)
            for i in range(r):
                for j in range(r):
                    Mn[i, j] = M[i, j] * mp.mpf(d - i - j)**n
            Rn = norm2(Si * Mn * Si.T)
            ratio_true = Rn / (t_true / theta)**n
            ratio_claim = Rn / wc**n
            worst = max(worst, ratio_true)
            assert ratio_true <= 1 + mp.mpf(10)**(-30), (trial, n, ratio_true)
        rows.append((float(al), float(theta), r, float(t_true), float(t_claim)))

    for row in rows:
        print('alpha=%.3f theta=%.2e r=%2d  t_true=%.4e  t_claim=4sqrt(r theta)=%.4e' % row)
    print('max over tests of R(s^n)/(t_true/theta)^n =', mp.nstr(worst, 6), ' (must be <= 1)')
    print('Gate B response bound R(s^n) <= w^n verified numerically (n=0..8).')
