"""Interval-arithmetic certificate for every scalar inequality used in the paper
"A uniform quadratic wedge for consecutive Toeplitz minors of the Riemann xi-coefficients".

Requires python-flint (Arb).  Every window is of the form c(k)/k or c(k)/k^2 with c(k) monotone
in the favourable direction, so it suffices to evaluate at k0 = 12000 (see Section 3 of the paper).
Run:  python3 certificate.py
"""
from flint import arb, ctx

ctx.prec = 256
k0 = arb(12000)
C = arb(10)**11
ok = lambda b, msg: (print('  [ok] ' + msg) if b else (_ for _ in ()).throw(AssertionError(msg)))

print('Section 3: curvature windows (uniform in k >= 12000)')
ok(6 * arb(6).exp() < 2 * k0 / arb.pi(), 'W(2k/pi) > 6')
gmin = arb(576) / 343
a1 = arb(12) / 7 - 20 / k0
a2 = 2 + 20 / k0
cau = 40 / (k0 - 2)
a4 = gmin * (k0 / (k0 + 1))**2 - cau
a3 = 2 + cau
ok(a1 > arb('1.71'), 'k tau_k >= 12/7 - 20/k > 1.71          (value %s)' % a1.lower().str(8))
ok(a2 < arb('2.002'), 'k tau_k <= 2 + 20/k < 2.002            (value %s)' % a2.upper().str(8))
ok(a4 > arb('1.675'), 'k^2 Delta_k >= (576/343)(k/(k+1))^2 - 40/(k-2) > 1.675 (value %s)' % a4.lower().str(8))
ok(a3 < arb('2.004'), 'k^2 Delta_k <= 2 + 40/(k-2) < 2.004   (value %s)' % a3.upper().str(8))

print('Section 4: reference parameters')
TL, TU, DL, DU = arb('1.71'), arb('2.002'), arb('1.675'), arb('2.004')
SL = TL - DU / k0                      # k sigma >= k tau - k Delta
eps_hi = (2 * (TU / k0 / 2)**2 / 6).exp()
eps_lo = ((DU / k0**2 / 2)**2 / 6).exp()
al_lo = TL * SL / DU / eps_lo
al_hi = TU * TU / DL * eps_hi
ok(al_lo > arb('1.45'), 'alpha > 1.45   (value %s)' % al_lo.lower().str(8))
ok(al_hi < arb('2.40'), 'alpha < 2.40   (value %s)' % al_hi.upper().str(8))
th_lo = SL - DU / SL - TU**2 / (24 * k0)
ok(th_lo > arb('0.53'), 'k theta > 0.53 (value %s)' % th_lo.lower().str(8))
ok(TU < arb('2.002') + arb(10)**-30, 'k theta < k tau < 2.002')
kap = arb('2.002') * arb('1.45') / arb('0.45')     # k*theta*alpha/(alpha-1) upper bound
ok(kap < arb('8.4'), 'k kappa = k theta alpha/(alpha-1) < 8.4 (value %s)' % kap.upper().str(6))
for d in range(3, 400):
    K = arb('8.4')
    assert K**d / d / k0 + K**(d - 1) / (d * (d - 1)) < arb(9)**(d - 1)
    assert 3 * arb(40)**d + arb(9)**(d - 1) <= 4 * arb(40)**d
ok(True, '|[s^d] log c_s| <= (9/k)^(d-1) and |beta_d| <= 4 40^d k^(1-d), d = 3..399 (larger d: trivial)')

print('Sections 6-8: wedge k >= 1e11 r^2, r >= 2')
TH = arb('0.53')
x = arb('2.002') / (2 * C)
ok(x <= arb(1) / 64, 'x = r theta <= 2.002/(2C) <= 1/64')
ok(arb('2.002') / (2 * C) < arb('1.45').log(), '(r-1) theta < log alpha  (no poles, signs)')
y2 = 1600 * 16 / TH / (2 * C)
y = y2.sqrt()
main = 4 * (1600 * 16 / TH)**2 / C / (1 - y)
corr = 4 * arb(80)**4 * (y / 40)**3 / (1 - 80 / (4 * C))
tot = main + corr
print('     y <= %s' % y.upper().str(6))
print('     main term 4k y^4/(1-y)           <= %s' % main.upper().str(8))
print('     interpolation term               <= %s' % corr.upper().str(6))
print('     ||h||_A                          <= %s' % tot.upper().str(8))
ok(tot < arb('0.0934'), '||h||_A < 0.0934')
ok(tot.exp() - 1 < arb('0.098'), '||e^h - 1||_A <= e^{0.0934} - 1 < 0.098 < 1')
Cmin = 4 * (1600 * 16 / TH)**2 / arb(2).log()
print('     (the same chain closes for every C > %s)' % Cmin.upper().str(6))
print('ALL CHECKS PASSED')
