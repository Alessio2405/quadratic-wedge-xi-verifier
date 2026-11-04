"""Plausibility check of Gate A and of the reference-coefficient bound, by Cauchy integrals.

f(z) = log I(z) - log Gamma(2z+1),  I(z) = int_0^oo u^{2z} Phi(u) du   (Michalowski's analytic log a).
We extract Taylor coefficients of f at k on |s| = rho, and of log c_s (analytic continuation of the
reference) on the same circle; then beta_d = [s^d](f - log c) for d>=2 (d>=4 matter).
"""
import mpmath as mp

mp.mp.dps = 40


def Phi(u):
    return sum((2 * mp.pi**2 * n**4 * mp.e**(9 * u) - 3 * mp.pi * n**2 * mp.e**(5 * u))
               * mp.e**(-mp.pi * n**2 * mp.e**(4 * u)) for n in range(1, 3))


def logI(z, u0):
    # integrand scaled by its value near the real saddle to avoid overflow; return log I
    lscale = 2 * z * mp.log(u0) + mp.log(Phi(u0))
    g = lambda u: mp.e**(2 * z * mp.log(u) + mp.log(Phi(u)) - lscale)
    pts = [u0 - 0.6, u0 - 0.3, u0 - 0.15, u0 - 0.05, u0, u0 + 0.05, u0 + 0.15, u0 + 0.3, u0 + 0.8]
    val = mp.quad(g, pts, maxdegree=10)
    return mp.log(val) + lscale


def saddle(k):
    # maximize 2k log u + log Phi(u)
    return mp.findroot(lambda u: 2 * k / u + mp.diff(lambda x: mp.log(Phi(x)), u), 1.0 + 0.25 * mp.log(k))


def coeffs(F, rho, N, dmax, unwrap=True):
    vals = [F(rho * mp.expjpi(2 * mp.mpf(j) / N)) for j in range(N)]
    # unwrap imaginary parts (continuous log along the circle)
    for j in (range(1, N) if unwrap else []):
        while vals[j].imag - vals[j - 1].imag > mp.pi:
            vals[j] -= 2j * mp.pi
        while vals[j].imag - vals[j - 1].imag < -mp.pi:
            vals[j] += 2j * mp.pi
    assert (not unwrap) or abs(vals[0].imag - vals[-1].imag) < 1, 'winding'
    out = []
    for d in range(dmax + 1):
        acc = mp.mpc(0)
        for j in range(N):
            acc += vals[j] * mp.expjpi(-2 * mp.mpf(j) * d / N)
        out.append((acc / N / rho**d))
    return out


def qp_int(a, qq, n):
    if n >= 0:
        p = mp.mpf(1)
        for j in range(n):
            p *= 1 - a * qq**j
        return p
    p = mp.mpf(1)
    for j in range(1, -n + 1):
        p *= 1 - a * qq**(-j)
    return 1 / p


for k in [12000, 100000]:
    u0 = saddle(k)
    f = lambda s: logI(k + s, u0) - mp.loggamma(2 * (k + s) + 1)
    # integer values -> u, v, alpha, theta
    fk = {s: f(mp.mpf(s)).real for s in (-1, 0, 1, 2)}
    tau = -(fk[-1] - 2 * fk[0] + fk[1]); sig = -(fk[0] - 2 * fk[1] + fk[2])
    u, v = mp.e**-tau, mp.e**-sig
    qq = (1 - v) / (v * (1 - u)); al = (1 - u) * (1 - v) / (v - u); th = mp.log(qq)
    z = 1 / al; Qr = 1 / qq
    # analytic log c_s
    def logc(s):
        ser = mp.fsum(z**m * (1 - mp.e**(-m * th * s)) / (m * (1 - mp.e**(-m * th))) for m in range(1, 400))
        return s * mp.log(1 - z) - th * s * (s - 1) / 2 + ser
    for s in range(-4, 5):
        assert abs(mp.e**logc(mp.mpf(s)) / ((1 - al)**s / qp_int(al, qq, s)) - 1) < mp.mpf(10)**-25
    rho = mp.mpf(k) / 50
    N = 64
    # f(k+s): subtract linear part analytically irrelevant for d>=2
    f0 = fk[0]
    f1 = (fk[1] - fk[-1]) / 2
    f2 = (fk[1] - 2 * fk[0] + fk[-1]) / 2
    # subtract an exact polynomial (any quadratic) to keep phases small; d>=3 coefficients unaffected
    F = lambda s: f(s) - (f0 + f1 * s + f2 * s**2)
    fc = coeffs(F, rho, N, 9)
    fc[2] += f2; fc[1] += f1; fc[0] += f0
    cc = coeffs(logc, rho, N, 9, unwrap=False)
    print('k=%d: k*tau=%.5f  alpha=%.5f  k*theta=%.5f ;  analytic log c_s matches integer c_s on [-4,4]' % (k, k * tau, al, k * th))
    print('   d   |f_d|*k^(d-1)   GateA 3*40^d   |[s^d]log c|*k^(d-1)  (9)^(d-1)   |beta_d|*k^(d-1)')
    for d in range(2, 10):
        fd = abs(fc[d].real) * mp.mpf(k)**(d - 1)
        cd = abs(cc[d].real) * mp.mpf(k)**(d - 1)
        bd = abs((fc[d] - cc[d]).real) * mp.mpf(k)**(d - 1)
        print('  %d   %.6e    %.2e     %.6e         %.2e   %.6e' % (d, fd, 3 * 40**d, cd, 9**(d - 1), bd))
