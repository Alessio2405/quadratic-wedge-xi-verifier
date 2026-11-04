"""How does R(s^n) and R(P_4) actually scale with r and theta for the q-Pochhammer reference?
Compare with the proven majorant (t/theta)^n, t = 4 sqrt(r theta)."""
import mpmath as mp
from gateB_numeric import qp, c, ldl, norm2

mp.mp.dps = 200

def responses(al, theta, r, polys):
    qq = mp.e**theta; d = r - 1
    M = mp.matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = c(d - i - j, al, qq)
    L, D = ldl(M)
    S = L * mp.diag([mp.sqrt(abs(x)) for x in D]); Si = S**-1
    out = []
    for p in polys:
        Mn = mp.matrix(r, r)
        for i in range(r):
            for j in range(r):
                Mn[i, j] = M[i, j] * p(mp.mpf(d - i - j))
        out.append(norm2(Si * Mn * Si.T))
    return out

P4 = lambda s: s * (s + 1) * (s - 1) * (s - 2)
polys = [lambda s: s, lambda s: s**2, lambda s: s**3, lambda s: s**4, P4]
al = mp.mpf('1.7777')
for k in [3000, 12000, 48000]:
    theta = mp.mpf('0.7655') / k
    print('k=%d' % k)
    print('   r    R(s)/sqrt(r/th)  R(s^2)/(r/th)  R(s^3)/(r/th)^1.5  R(s^4)/(r/th)^2  R(P4)/(r/th)^2   majorant 16^2')
    for r in [4, 8, 16, 32]:
        R = responses(al, theta, r, polys)
        z = r / theta
        print('%4d    %.4f          %.4f         %.4f             %.4f          %.4f' % (
            r, R[0] / mp.sqrt(z), R[1] / z, R[2] / z**1.5, R[3] / z**2, R[4] / z**2))
