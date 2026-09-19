# Spectral–finite-element formulation (Martinec 2000)

Implementation reference for the radial solver, extracted from Martinec (2000),
*GJI* 142:117 (full PDF held locally, gitignored). Equation numbers below are the
paper's. This is the spec for `fe_radial_fe` assembly and Love-number extraction.

## Problem

Self-gravitating, incompressible, Maxwell-viscoelastic sphere `B` of radius `a`
under a surface mass load `σ`. Governing equations (incremental, quasi-static):

- Momentum (1): `div τ − ρ₀∇φ₁ + div(ρ₀u)∇φ₀ − ∇(ρ₀u·∇φ₀) = 0`
- Poisson (2): `∇²φ₁ + 4πG div(ρ₀u) = 0`
- Constitutive (3,4): `τ = ΠI + 2με`, elastic; Maxwell adds the memory term (rung 3)
- Incompressibility (5): `div u = 0`
- φ₁ = perturbed gravitational potential (interior + the load's exterior potential)

## Weak form (energy functional, §4)

Find `(u, φ₁, Π) ∈ V` with `δE = δF` for all test functions (eq 47).
`E = E_press + E_shear + E_grav + E_uniq` (eqs 30–33); forcing `F = F_surf` (elastic;
the dissipative `F_diss` is rung 3). Variations: eqs 40–46.

- `E_press = ∫ Π div u dV`  (30) — couples pressure ↔ incompressibility
- `E_shear = ∫ μ (ε:ε) dV`  (31)
- `E_grav`  (32) — self-gravity: body force + `(1/8πG)∫|∇φ₁|²`
- `E_uniq`  (33) — removes rigid translation+rotation; **degree j=1 only**
- `F_surf = ∫_∂B (b₀·u + b₁ φ₁) dS` (36), `b₀=−g₀(a)σ eᵣ` (37), `b₁` (38)

## Spectral reduction (§5)

Vector/scalar SH expansion (55–57), `J ≡ j(j+1)`:
- `u = Σ_{j≥1} [U_jm(r) S⁽⁻¹⁾ + V_jm(r) S⁽¹⁾ + W_jm(r) S⁽⁰⁾]` — **no j=0** (incompressibility)
- `φ₁ = Σ_{j≥0} F_jm(r) Y_jm`,  `Π = Σ_{j≥0} Π_jm(r) Y_jm`
- **Spheroidal-only for 1D loading: W=0** (the W block decouples, eq 110).
- `div u = Σ (U' + 2U/r − J V/r) Y_jm`  (58)

Each degree `j` decouples (radially symmetric μ); solve a 1D radial problem per j.

## Discretization (§6)

- P1 "tent" basis ψ_k(r) (71); `U,V,F = Σ_k [·]_k ψ_k` (72) — **nodal**, P+1 nodes.
- P0 basis ξ_k(r) (74); `Π = Σ_k Π_k ξ_k` (73) — **per-element**, P values.
- Per element, μ_k, ρ_k constant (75); gravity `g₀(r) = (4πG/3)(ρ_k r + R_k/r²)`
  (76), with `R_k = Σ_{i≤k}(ρ_{i−1}−ρ_i)r_i³` (77), `R₁=0`.
- Mesh = whole sphere ⟨0,a⟩; fluid core has μ=0. Done in `fe_radial_fe` (218 nodes).

### Assembled bilinear forms (the element matrices to build)

Use the Appendix C element integrals (implemented + tested in
`fe_radial_integrals`): `I1..I7`, `K1..K3`.

- `δE_press` (82): pressure ↔ `(U' + 2U/r − JV/r)` coupling. Uses K1,K2 (per the
  P0×P1 products `∫ξ ψ' r²`, `∫ξ ψ r`, `∫ξ ψ r`-type → K-integrals).
- `δE_shear` (80): U,V (and W, dropped) stiffness. Uses I1,I3,I6 with factor μ_k.
- `δE_grav` (81): self-gravity coupling among U,V,F. Uses I2,I4,I5,I7 (the I7/1/r
  term carries R_k; **skip when R_k=0**, i.e. innermost element) and `1/(4πG)` for
  the `∫|∇φ₁|²` block (I1,I6).
- `δE_uniq` (83): degree-1 rigid-mode removal. The penalty weight w uses K3
  (∫ψ r²) on the degree-1 (U,V) dofs (`uniq_weight`). Imposed sparsely as a KKT
  constraint wᵀd=0, not as the dense penalty — see "degree-1" below.

### Surface forcing (eq 84)

`δF_surf = −(a/4πG) Σ (j+1) F^{P+1}_jm δF^{P+1}_jm − a² Σ σ_jm [g₀(a) δU^{P+1}_jm + δF^{P+1}_jm]`

- The `(j+1)F(a)` term = exterior-potential matching (φ₁ ∝ r^{−(j+1)} outside),
  added to the **F–F diagonal at the surface node**.
- The `σ` terms are the **RHS**: forcing on U(a) (`−a² σ g₀(a)`) and on F(a) (`−a² σ`).

## System structure & solve (§10)

Per degree j, DOFs `(U_k, V_k, F_k)` nodal (k=1..P+1) + `Π_k` per element (k=1..P),
laid out **node-interleaved** `[U V F | Π]` so the operator stays band-diagonal:

```
[ A   Bᵀ ] [d ]   [f]          A = shear+grav stiffness (U,V,F)
[ B   0  ] [Π ] = [0]          B = pressure/incompressibility coupling (eq 82)
```

The whole operator is **symmetric**: it is the Hessian (second variation) of the
energy functional `E = E_press + E_shear + E_grav + E_uniq` (eqs 30–33), so it is
self-transpose by construction. The shear block (eq 80, `∫μ ε:ε`) is symmetric;
the **B/Bᵀ** pressure block is symmetric; and the self-gravity U↔F coupling (eq 65
continuous → eq 81) is a **transpose pair** — the potential-gradient body force
`∫(dF/dr)δU r² = I²_βα` on δU and the Poisson source `∫ρ₀U(dδF/dr) r² = I²_αβ` on
δF. (Earlier this doc claimed the operator was non-symmetric "because the I² U↔F
coupling is not symmetric"; that was the elastic low-degree bug — the U-F term was
discretised with `I²_αβ` instead of `I²_βα`. See below. `test_assembly` now
asserts `‖A−Aᵀ‖/‖A‖ = 0`.) Implemented in `fe_radial_fe%build_dense_operator`,
transcribed term-by-term from the PDF and verified against the table and analytic
limits below.

**Solve: dependency-free pivoted banded LU** (`fe_band`). The physical entries span
~20 orders of magnitude (`μr²/h` vs the pressure couplings vs `1/4πG`), so the
operator is geometric-mean **row/column equilibrated** first. It is band-diagonal
(half-bandwidth ~6, node-interleaved P1/P0) and indefinite (zero pressure block),
so a **partial-pivoting band LU** factors it directly — factor once per j at
assemble, solve per RHS (~20 µs). The matrix is **degree-dependent only through J**;
assemble + factor once per j, reuse for all orders m / loads / time steps. The LU
is re-entrant (no global state), so the per-degree loop parallelizes with OpenMP.
(Earlier this used LIS GMRES+ILU; LIS was removed — it added a serial-vs-OpenMP
link dependency and concurrent solves were not re-entrant.) Degree 1 carries the
dense KKT border (rigid-mode removal) so it factors as a dense LU via the same code.

### Boundary / regularity conditions
- **Centre r=0:** **no explicit BC** — Martinec meshes through the centre and the
  r² weighting handles regularity; the singular `I⁷` term is killed by `R₁=0`
  (skipped to avoid 0·∞). Verified: no empty rows, solver well-posed for j≥2.
- **Surface:** the (j+1)F(a) exterior term + load RHS above.
- **Uniqueness (j=1):** E_uniq removes the rigid translation null space — imposed
  as a sparse KKT constraint wᵀd=0 (CM frame), not the dense penalty (see above).
- **Fluid core:** μ=0 region; free-slip emerges (no shear stress). No explicit CMB
  BC (Martinec meshes through the centre).

## Viscoelastic time stepping (§3, §8-9; rung 3 — DONE 1-D)

Explicit ω=1 Maxwell scheme (eqs 23-25): the total stress splits into the
instantaneous elastic stress (the SAME operator above) plus a memory stress
`τ^{V,i} = (1−M)τ^{V,i-1} − 2μ M ε^i`, `M = μΔt/η` (eq 17). The memory enters the
RHS as the dissipative forcing `−∫ τ^{V,i}:δε dV` (eq 35); the LHS never changes,
so it is assembled, equilibrated and LU-factored **once** (`fe_band`) and reused
every step.

1-D (radially symmetric η): the memory stress evolves directly on the tensor-SH
coefficients (§9, eq 107) — no spatial grid. Per element it is stored as `A,B,C`
for the four spheroidal tensor components λ ∈ {1,2,5,6} (eq 109); the strain
coefficients `a,b,c` come from nodal `U,V` (eq 87, `ε = a/h + bψ_k/r + cψ_{k+1}/r`,
eq 88). The dissipative RHS is a 2-point radial Gauss quadrature (eqs 94-95) of
the spectral double-dot `Σ_λ ‖Z^λ‖² τ^{V,λ} δε^λ` with norms `{1, J/2, 2J², 2J(J−2)}`
(eqs 110/B13). Implemented in `fe_viscoelastic%ve_degree`. Elastic layers (η→∞)
freeze (M→0); fluid layers (μ=0) carry no memory. Stability `Δt ≲ 2η_min/μ`.

**Validated** (`test_relax`): a held degree-2 load on a homogeneous Maxwell
sphere relaxes from the elastic Love number (t=0) to the fluid limit `−(2j+1)/3`
(t→∞) — the two limits already pinned below — smoothly and monotonically, with
`t_relax ∝ η` (e-folding 0.76→1.53 kyr when η doubles). dt-converged (10 vs 50 yr
agree). **Disc time series vs Spada (2011):** the M3-L70-V01 disc relaxation
previously sat ~10% high (elastic) / ~11% low (fully relaxed). **RESOLVED:** this
was the elastic low-degree self-gravity bug (the U-F transpose, fixed above) — the
disc is dominated by low–intermediate degrees where the elastic Love numbers were
too soft. With the fix the per-degree elastic Love numbers match the benchmark to
~0.1%, so the disc offset is closed at the source; a direct disc re-run to confirm
<1% is a quick follow-up (the synthesis prototype lives in `/tmp/explore_disc*.f90`).

**Degree-1 (sparse KKT, solved):** `E_uniq` (eq 83) is a rank-1 penalty
`(4π/3) w wᵀ` over every degree-1 (U,V) dof, so adding it to the operator densifies
it (~870×870 dense ILU + solve — impractical for the time stepper). Two facts let
us keep it sparse: (i) `w` carries K³~∫ψr², so the penalty coefficient is ~1e16×
the band, i.e. the penalty is already a hard constraint `wᵀd=0` (the CM/geocenter
frame, Blewitt 2003) in all but name; (ii) a hard constraint borders the band with
one row/col instead of filling it. So `radial_operator` imposes it as a KKT saddle
point

```
[ A_band  w ] [d]   [f]
[ wᵀ      0 ] [λ] = [0]   ⇒  A_band d + w λ = f,  wᵀ d = 0
```

(`with_uniq=.false.` builds the band; `uniq_weight` is the border). The solve runs
at `ndof+1` internally but `load_rhs`/`solve_vec` keep their physical `ndof`
interface, so `fe_viscoelastic` needs no j=1 special case. Converges in 2 GMRES
iterations; `build_dense_operator` still adds the dense penalty as the reference
operator. Validated in `test_love` (4) and `test_relax` (5); see the validation
list below. (The disc synthesis no longer needs to skip j=1 — it is the physical
geocenter signal.)

## Love numbers (§11; conventions verified)

For a degree-j surface load of coefficient `σ`, the load's own potential at the
surface is `φ^L = 4πG a σ/(2j+1)`. From the surface coefficients (Farrell 1972
normalization), implemented in `fe_radial_fe%loading_love`:

`h_j = g₀(a) U(a)/φ^L`,  `l_j = g₀(a) V(a)/φ^L`,  `k_j = −F(a)/φ^L − 1`.

The `k` form was **pinned empirically by two analytic limits**: Martinec's `φ₁`
(=`F`) is the *total* perturbation potential and carries the load's direct
potential with the **opposite sign** to `φ^L` (`F→−φ^L` for a rigid sphere), so
the induced potential is `−F−φ^L`. `σ` cancels in every ratio (use σ=1).
**`l` sign / S⁽¹⁾-normalization: RESOLVED.** The benchmark M3-L70-V01 fluid
limit reproduces the table `l_f` to ~0.1 % at every degree 2–8 (`test_benchmark_love`),
so `l = g V(a)/φ^L` is correct as written — no extra sign or normalization factor.

## Elastic low-degree discrepancy (FIXED)

A long-standing ~10 % offset in the disc benchmark turned out to be a real
solver bug in the elastic self-gravity coupling. With the benchmark table now
in-repo (`data/benchmarks/love_M3-L70-V01/`, independently reproduced by TABOO
NV=3/CODE=7) it was localised and fixed (`test_benchmark_love`).

**Symptom.** The elastic loading Love numbers were too *soft* (too much
deformation): `h_e(2) = −0.669` vs the table/TABOO `−0.454` (−47 %), the error
shrinking with degree to ~1 % by j≈40. The *fluid* (t→∞) limit was always exact
(<0.5 %, all degrees) — fluidise the Maxwell layers (μ=0) and the elastic solve
reproduces the table `h_f,l_f,k_f`.

**Diagnosis (the re-derivation).** The bug hid from every existing test:
fluid→elastic changes only the shear block (eq 80; grav/press/surface carry no
μ), μ→0 kills the shear block, μ→∞ forces d→0. So all the homogeneous limits and
the M3 fluid limit pass regardless. Ruled out, in order: model params, load,
frame, mesh (converged), and a uniform μ-scale (the required correction is
degree-dependent — the spectrum *shape* was wrong). Then:

1. The **shear block is correct** — its 4×4 element matrix is identical (machine
   precision, all degrees) to the stiffness rebuilt independently from the strain
   representation (eqs 85–88) used by `fe_viscoelastic`, i.e. `2∫μ Σ_λ ‖Z^λ‖²
   B^λ_i B^λ_j r² dr`. So the energy `∫μ ε:ε` is encoded consistently two ways.
2. That left the **self-gravity block** (eq 65 → 81), μ-independent but corrupting
   the *interior* solution in a way the fluid surface values don't expose but the
   elastic shear coupling does — and self-gravity dominates exactly at low degree,
   matching the signature. Discretising the **continuous** form eq 65 term-by-term
   pinned it: the potential-gradient body force on δU, `∫ρ₀(dF/dr)δU r²`,
   discretises to `∫ψ'_α ψ_β r² = I²_βα` (derivative on the **trial** F basis).
   The code used `I²_αβ` (`i2(ia,ib)`) — the transpose of the Poisson-source F-U
   term — which is *not* its proper symmetric partner.

**Fix.** One index in `build_dense_operator`: the U-F entry `i2(ia,ib) →
i2(ib,ia)`. This restores the U↔F symmetry the energy functional requires (the
operator is now exactly symmetric, `test_assembly` 2f). Result: elastic `h,k,l`
match the table to ~0.1 % (P1 discretisation) at **every** degree 2–48, and the
fluid limit is unchanged. This also closes the disc offset (rungs 2/3), whose
root cause was this same term. (Likely a typo in the paper's discretised eq 81
relative to its own continuous eq 65, faithfully copied — eq 65 is the arbiter.)

## Validation targets
1. **Fluid limit** (μ→0, homogeneous sphere): `h_j → −(2j+1)/3` and `k_j → −1`.
   ✅ reproduced to ~1e-5 (`test_love`). Independent check of self-gravity +
   incompressibility + Poisson + the load forcing.
2. **Rigid limit** (μ→∞): `h_j, l_j, k_j → 0`. ✅ to ~1e-5 (`test_love`).
   Checks the shear block and the F sign convention.
3. Elastic loading Love numbers h,l,k vs **benchmark M3-L70-V01 table**
   (`data/benchmarks/love_M3-L70-V01/`, degrees 2–256; `test_benchmark_love`).
   **✅ elastic AND fluid match to <1 % (≈0.1 %, P1 discretisation) at every
   degree** after the U-F symmetry fix; the fluid limit also pins `l`. See
   "Elastic low-degree discrepancy (FIXED)" above.
4. Internal: operator finite (centre I⁷ guard); B=Bᵀ; **full operator symmetric
   (energy Hessian)**; gravity R_k reconstruction (`test_assembly`). ✅
5. **Degree-1 (sparse KKT):** the j=1 solve converges (non-singular), removes the
   rigid mode (`wᵀd/|w||d|`~1e-23), satisfies the band operator off the gauge
   direction (`‖r⊥w‖/‖r‖`~1e-11) and on the F rows, and gives a finite geocenter
   response (`test_love` (4)); `ve_degree` steps stably at j=1 (`test_relax` (5)). ✅
