"""Independent symbolic checks (sympy, generic symbols) of the reference algebra."""
import sympy as sp

u, v, al, q, B, Q = sp.symbols('u v alpha q B Q', positive=True)


def qp(a, qq, n):
    """(a;qq)_n for integer n with the standard extension (a;q)_{-m} = 1/(a q^{-m}; q)_m."""
    if n >= 0:
        return sp.prod([1 - a * qq**j for j in range(n)]) if n > 0 else sp.Integer(1)
    m = -n
    return 1 / sp.prod([1 - a * qq**(-j) for j in range(1, m + 1)])


def c(s, a=al, qq=q):
    return (1 - a)**s / qp(a, qq, s)

# ---- 1. four-point matching, with alpha, q as functions of u, v
qs = (1 - v) / (v * (1 - u))
als = (1 - u) * (1 - v) / (v - u)
for s, target in [(-1, u), (0, 1), (1, 1), (2, v)]:
    assert sp.simplify(c(s).subs({al: als, q: qs}) - target) == 0, s
print('[1] matching c_{-1}=u, c_0=c_1=1, c_2=v : OK (symbolic in u,v)')

# sinh forms of alpha and theta
tau, sig = sp.symbols('tau sigma', positive=True)
al_sinh = 2 * sp.sinh(tau / 2) * sp.sinh(sig / 2) / sp.sinh((tau - sig) / 2)
q_exp = sp.exp((tau + sig) / 2) * sp.sinh(sig / 2) / sp.sinh(tau / 2)
sub = {u: sp.exp(-tau), v: sp.exp(-sig)}
for expr, ref in [(als.subs(sub), al_sinh), (qs.subs(sub), q_exp)]:
    for tv, sv in [(sp.Rational(3, 1000), sp.Rational(29, 10000)), (sp.Rational(1, 7), sp.Rational(1, 9))]:
        val = (expr - ref).subs({tau: tv, sig: sv}).evalf(60)
        assert abs(val) < 1e-50
print('[1b] sinh representations of alpha and q: OK (60-digit spot checks)')

# ---- 2. cocycle with negative indices, (a;q)_{m+n} = (a;q)_m (a q^m; q)_n for all integers m,n
for m in range(-5, 6):
    for n in range(-5, 6):
        lhs = qp(al, q, m + n)
        rhs = qp(al, q, m) * qp(al * q**m, q, n)
        assert sp.simplify(sp.together(lhs - rhs)) == 0, (m, n)
print('[2] integer-index cocycle (a;q)_{m+n}=(a;q)_m(aq^m;q)_n for m,n in [-5,5]: OK')

# ---- 3. congruence M_ij = c_{d-i-j} = c_d A_i U_ij A_j, with Q=1/q, B=alpha q^{d-1}
for r in range(2, 7):
    d = r - 1
    Bv = al * q**(d - 1)
    Qv = 1 / q
    for i in range(r):
        for j in range(r):
            n = i + j
            lhs = qp(al, q, d) / qp(al, q, d - n)
            rhs = qp(Bv, Qv, n)
            assert sp.simplify(sp.together(lhs - rhs)) == 0
            Ai = qp(Bv, Qv, i) / (1 - al)**i
            Aj = qp(Bv, Qv, j) / (1 - al)**j
            U = qp(Bv, Qv, i + j) / (qp(Bv, Qv, i) * qp(Bv, Qv, j))
            assert sp.simplify(sp.together(c(d - n) - c(d) * Ai * U * Aj)) == 0
print('[3] reversed-Hankel congruence incl. negative indices, r=2..6: OK (symbolic alpha,q)')

# ---- 4. LDL with generic B, Q
def qbin(n, k, QQ):
    if k < 0 or k > n:
        return 0
    return sp.cancel(qp(QQ, QQ, n) / (qp(QQ, QQ, k) * qp(QQ, QQ, n - k)))

N = 7
for i in range(N):
    for j in range(N):
        U = qp(B, Q, i + j) / (qp(B, Q, i) * qp(B, Q, j))
        S = sum(qbin(i, m, Q) * qbin(j, m, Q) * B**m * Q**(m * (m - 1)) * qp(Q, Q, m) / qp(B, Q, m)
                for m in range(min(i, j) + 1))
        assert sp.cancel(sp.together(U - S)) == 0, (i, j)
print('[4] U = L D L^T with L=[i,j]_Q, D_m = B^m Q^{m(m-1)} (Q;Q)_m/(B;Q)_m : OK for i,j<7 (generic B,Q)')

# ---- 5. R_1 = L^{-1} diag(Q^i) L is lower bidiagonal with stated entries
N = 7
L = sp.Matrix(N, N, lambda i, j: qbin(i, j, Q))
R1 = sp.simplify(L.inv() * sp.diag(*[Q**i for i in range(N)]) * L)
for i in range(N):
    for j in range(N):
        if i == j:
            exp_ = Q**i
        elif i == j + 1:
            exp_ = -Q**(i - 1) * (1 - Q**i)
        else:
            exp_ = 0
        assert sp.simplify(R1[i, j] - exp_) == 0
print('[5] R_1 lower bidiagonal: (R1)_mm=Q^m, (R1)_{m,m-1}=-Q^{m-1}(1-Q^m): OK (N=7)')

# ---- 6. pivot ratio and whitened subdiagonal
m = sp.symbols('m', integer=True, positive=True)
Dm = lambda mm: B**mm * Q**(mm * (mm - 1)) * qp(Q, Q, mm) / qp(B, Q, mm)
for mm in range(1, 7):
    ratio = sp.simplify(Dm(mm) / Dm(mm - 1))
    expected = -B * Q**(2 * mm - 2) * (1 - Q**mm) / (B * Q**(mm - 1) - 1)
    assert sp.simplify(ratio - expected) == 0
    # whitened entry squared: (R1)_{m,m-1}^2 * |D_{m-1}|/|D_m|
    w2 = (Q**(mm - 1) * (1 - Q**mm))**2 / (B * Q**(2 * mm - 2) * (1 - Q**mm) / (B * Q**(mm - 1) - 1))
    assert sp.simplify(w2 - (1 - Q**mm) * (B * Q**(mm - 1) - 1) / B) == 0
print('[6] pivot ratio D_m/D_{m-1} = -BQ^{2m-2}(1-Q^m)/(BQ^{m-1}-1) and whitened subdiagonal^2 =(1-Q^m)(BQ^{m-1}-1)/B: OK')

# ---- 7. Lagrange basis on nodes -1,0,1,2
s = sp.symbols('s')
nodes = [-1, 0, 1, 2]
def ell(k):
    return sp.expand(sp.prod([(s - t) / (k - t) for t in nodes if t != k]))
assert ell(-1) == sp.expand(-s**3/6 + s**2/2 - s/3)
assert ell(1) == sp.expand(-s**3/2 + s**2/2 + s)
assert ell(2) == sp.expand(s**3/6 - s/6)
w = sp.symbols('w', positive=True)
worst = 0
for dd in range(4, 40):
    I3 = sp.expand(sum(sp.Integer(t)**dd * ell(t) for t in nodes))
    norm = sum(abs(I3.coeff(s, e)) * w**e for e in range(4))
    # check norm <= 2^d w^3 for all w>=1 : norm/w^3 is decreasing in w, so check at w=1
    ratio = norm.subs(w, 1) / 2**dd
    worst = max(worst, ratio)
    assert ratio <= 1
print('[7] Lagrange basis OK;  max_d ||I_3[s^d]||_A /(2^d w^3) at w=1 over d=4..39 =', sp.nsimplify(worst), '(<=1): OK')
