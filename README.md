# Quadratic wedge for Toeplitz minors of the Riemann ξ-coefficients — verifier

Companion code for the paper

> **Alessio Doria – Independent Researcher** (numero2405@gmail.com), *A uniform quadratic wedge for consecutive Toeplitz minors of the Riemann ξ-coefficients* (2026).
> PDF and LaTeX source in [`paper/`](paper/).

**Main Theorem.** Let ⅛ ξ(½ + ix/2) = Σ (−1)^k a_k x^{2k} and D_{r,k} = det[a_{k+j−i}]_{i,j=0}^{r−1}.
Then D_{r,k} > 0 for every r ≥ 2 and every k ≥ 10¹¹ r².

This improves the cubic wedge k ≥ 10¹⁸ r³ of W. Michałowski ([arXiv:2607.16795](https://arxiv.org/abs/2607.16795))
by replacing the three-point log-quadratic comparison sequence with a four-point q-Pochhammer reference
c_s = (1−α)^s / (α; 𝔮)_s, whose reversed Hankel block has an exact q-Pascal LDLᵀ factorization.

## What this repository verifies

| Script | Role in the proof | What it checks |
|---|---|---|
| `verifier/certificate.py` | **Required** | Every scalar inequality of §§3–8 in Arb ball arithmetic (256 bits), uniformly in k ≥ 12000 and in the wedge k ≥ 10¹¹ r²: curvature windows, α and θ windows, reference-coefficient bound, ‖h‖_𝒜 < 0.0934, ‖e^h − 1‖_𝒜 < 0.098. |
| `verifier/sym_check.py` | Regression | Symbolic (generic parameters) checks of the four-point matching, integer-index q-Pochhammer cocycle, reversed-Hankel congruence (incl. negative indices), q-Pascal LDLᵀ, bidiagonal dilation, whitened sub-diagonal, Lagrange interpolation norms. |
| `verifier/gateB_numeric.py` | Regression | The all-degree response bound R(sⁿ) ≤ (t/θ)ⁿ, with S obtained from a *numerical* LDLᵀ at 260 digits (independent of the closed-form pivots). |
| `verifier/response_scaling.py` | Illustration | Actual growth of R(sⁿ) in the variable r/θ. |
| `verifier/taylor_check.py` | Illustration | Taylor coefficients of log a and of the reference at k = 12000, 10⁵ via Cauchy integrals; size of β_n. |
| `verifier/xi_data.py` | Illustration | τ_k, Δ_k, α, θ on the true ξ-coefficients. |
| `verifier/whitened_experiment.py` | Illustration | Actual whitened perturbation, three-point vs four-point reference. |

The proof itself is in the paper. The symbolic and numerical scripts are regression tests and illustrations; only
`certificate.py` is part of the logical chain. The two analytic inputs — Michałowski's all-degree Taylor bound
(Gate A, [arXiv:2607.16795](https://arxiv.org/abs/2607.16795), Prop. 3.7) and the Liu–Mao log-curvature estimate
([arXiv:2610.04972](https://arxiv.org/abs/2610.04972), Lemma 3.3) — are imported, not re-certified here.

## How to run

```bash
pip install -r requirements.txt
./run_all.sh               # certificate + symbolic + numerical checks (~2 min)
./run_all.sh --with-data   # also downloads the Liu–Mao certified λ_n data and runs §9
```

Outputs go to `output/`; reference outputs are in `expected_output/`.
`--with-data` clones [YXLiu97/riemann-xi-infinite-log-concavity](https://github.com/YXLiu97/riemann-xi-infinite-log-concavity)
into `data/` (not redistributed here).

## Use of AI tools

The following AI models were used in the preparation of this work: Astra (GPT, OpenAI) and Claude (Anthropic).

## License

Code: MIT (see `LICENSE`). Paper: © the author.
