"""Load Liu–Mao certified lambda intervals, cross-check independently, and compute tau, Delta, alpha, theta."""
import json, sys
import mpmath as mp

mp.mp.dps = 110
PATH = sys.argv[1]


def dec(e):
    sgn, man, ex, _ = e
    return (-1)**sgn * mp.mpf(man) * mp.mpf(2)**ex


def load():
    d = json.load(open(PATH))
    lam = []
    width = mp.mpf(0)
    for lo, hi in d['lambda_intervals']:
        a, b = dec(lo), dec(hi)
        lam.append((a + b) / 2)
        if a > 0:
            width = max(width, (b - a) / a)
    return lam, width


if __name__ == '__main__':
    lam, width = load()
    print('n_max =', len(lam) - 1, ' max relative interval width =', mp.nstr(width, 3))

    # --- independent check of lambda_n via the moment integral (Polya kernel)
    def Phi(u):
        return sum((2 * mp.pi**2 * n**4 * mp.e**(9 * u) - 3 * mp.pi * n**2 * mp.e**(5 * u)) * mp.e**(-mp.pi * n**2 * mp.e**(4 * u)) for n in range(1, 4))

    mp.mp.dps = 60
    for n in [10, 1000, 12000]:
        m2n = mp.quad(lambda u: u**(2 * n) * Phi(u), [0, 0.5, 1, 1.5, 1.75, 2, 2.5, 4])
        lam_ind = 8 * mp.mpf(4)**n * m2n / mp.factorial(2 * n)
        print('n=%5d  lambda(LM data) / lambda(own integral) - 1 = %s' % (n, mp.nstr(lam[n] / lam_ind - 1, 5)))
    mp.mp.dps = 110

    def tau(k):
        return mp.log(lam[k]**2 / (lam[k - 1] * lam[k + 1]))

    def Qf(x):
        w = mp.lambertw(2 * x / mp.pi).real
        return 2 * w / (x * (1 + w)), w

    print('\n  k     k*tau_k     k*Q(k)     k^2|tau-Q|  (<=20?)   k^2*Delta_k  (in (5/3,2.1)?)   alpha     k*theta')
    for k in [1000, 5000, 10000, 11000, 12000, 12010, 12014]:
        t = tau(k); s = tau(k + 1); Dl = t - s
        q, w = Qf(mp.mpf(k))
        al = 2 * mp.sinh(t / 2) * mp.sinh(s / 2) / mp.sinh(Dl / 2)
        th = (t + s) / 2 + mp.log(mp.sinh(s / 2) / mp.sinh(t / 2))
        print('%6d  %.6f  %.6f  %10.4f   %10.5f   %10.5f  %10.5f' % (k, k * t, k * q, k * k * abs(t - q), k * k * Dl, al, k * th))
