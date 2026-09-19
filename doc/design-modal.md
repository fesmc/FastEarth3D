# FastEarth3D — modal viscoelastic response (design)

This document specifies a new viscoelastic response representation, `RESP_MODAL`,
that reduces the per-degree relaxation to a chosen number of normal modes `K`.
It is the design reference for the implementation on branch `modal`.

It complements `design.md`. The full Martinec time-domain FE method (`RESP_VE`)
is the *reference* this approximation is validated against and the limit it
converges to as `K → all`; it is **retained permanently** and remains selectable
(§7).

## 1. Motivation

FastEarth3D today has two endpoints with nothing between them:

- **Full VE (`RESP_VE`, Martinec FE)** — accurate, but pays a banded radial
  solve every step, carries a heavy radial memory-stress state
  (`Are/Aim/Bre/Bim/Cre/Cim(NLAM, ne, nk)`), is bounded by a Maxwell `Δt`
  stability ceiling (hence sub-stepping), and routes lateral viscosity through
  the expensive 6-component tensor-SH path (`fe_tensor_sh`).
- A hypothetical **single-mode** model (LV-ELVA, as in FastIsostasy) — cheap and
  unconditionally stable, but keeps only one relaxation mode: lossy, with no
  knob to recover accuracy.

The real Earth's loading response, per harmonic degree, is a **sum of a handful
of exponential relaxation modes**. The natural design is therefore neither
endpoint but a **dial**: choose how many modes `K` to keep per degree.

- `K = 1`   → the dominant mode ≈ LV-ELVA;
- `K = all` → **converges to** the full FE model (`RESP_VE`), to a set tolerance
  ("all" = all modes above a residue threshold, §2);
- `1 < K < all` → tune fidelity to the modelling context — fast continental
  deglaciation at low `K`, careful regional GIA at high `K` — in one code path,
  one validation harness, one coupling interface.

This is the **LV-ELVA representation translated into our global
spherical-harmonic framework**. We do *not* import FastIsostasy's regional
Cartesian / FFT machinery; only its physical idea (a depth-collapsed,
mode-limited relaxation with stacked lateral viscosity).

### Relation to the "no normal-mode root finding" constraint

`design.md §1` rules out normal-mode root finding — meaning the analytic
secular-determinant approach of SELEN/TABOO/ALMA. We do **not** do that. We
extract the relaxation modes from the **already-validated FE relaxation
propagator** (§5). This keeps the existing FE infrastructure as the single
source of truth: the modal basis is a property of the same discretized operator
the full model uses, so `K = all` converges to `RESP_VE`. It is a one-time init
cost, not a per-step root find.

## 2. Physical basis — how many modes

For an **incompressible**, self-gravitating, layered Maxwell Earth the
relaxation spectrum per degree is **finite and discrete** (this is exactly why
incompressibility is assumed — it keeps the count finite; compressibility would
make it a continuum). Roughly one mode per interface:

- **M0** — the fundamental mantle mode (dominant for surface loads);
- one **buoyancy (M-type)** mode per internal density discontinuity;
- one **transient (T)** mode per viscosity contrast / Maxwell layer;
- **C0** — the core–mantle-boundary mode.

For an M3-L70-V01-class structure (lithosphere + 3 mantle Maxwell layers + fluid
core) that is **~10 physical modes, of which typically 1–3 dominate** at any
given degree (which ones is degree-dependent). The FE memory space has ~500–600
DOFs per degree (§5), but the surplus modes are numerical with negligible
surface residue. Hence:

> **`n_modes = K` keeps the `K` highest-ranked modes; "all" keeps every mode
> with residue above a tolerance.**

Physical content saturates well before `K` reaches double digits.

**Mode ranking (`mode_rank` namelist):** the rank metric is a parameter, since
which is "best" is problem-dependent:

- `"isostatic"` (default) — `|r^u_k·b_k·τ_k|`: modes that dominate the final
  relaxed (fluid-limit) uplift; `K=1` ≈ M0;
- `"rate"` — `|r^u_k·b_k|`: weights by initial relaxation rate;
- `"residue"` — `|r^u_k|`: pure surface coupling.

All three are cheap byproducts of the solve, so `mode_rank` is sweepable.

## 3. Core idea

Per degree `l`, the FE-discretized viscoelastic system is linear time-invariant.
Its load→surface response, in modal form, is

```
u_lm(t) = u_el(l)·σ_lm  +  Σ_{k=1..K} r^u_k(l) · ξ_k,lm(t)
N_lm(t) = N_el(l)·σ_lm  +  Σ_{k=1..K} r^N_k(l) · ξ_k,lm(t)
V_lm(t) = V_el(l)·σ_lm  +  Σ_{k=1..K} r^V_k(l) · ξ_k,lm(t)

ξ̇_k,lm = −ξ_k,lm / τ_k(l)  +  b_k(l)·σ_lm     (one scalar ODE per mode per coeff)
```

- `u_el, N_el, V_el` — instantaneous **elastic** surface gains (the `t→0⁺`
  limit; same quantities `RESP_ELASTIC`/`RESP_VE` already compute).
- `{τ_k(l), r^{u,N,V}_k(l), b_k(l)}` — relaxation times, **modal residues**, and
  load projections from the per-degree solve (§5).
- The **isostatic (fluid) limit** is reproduced automatically:
  `u_∞ = u_el + Σ_k r^u_k b_k τ_k` etc. — a check, not a separately-imposed datum.

**State shrinks** from the full radial memory tensor to **`K` complex scalars per
(l,m)** — `ξ_k,lm`. Each obeys a scalar linear ODE.

**Exact time integration.** Over a step with constant load `σ`, each mode
advances by its analytic solution:

```
ξ_k ← e^{−Δt/τ_k} · ξ_k  +  b_k τ_k (1 − e^{−Δt/τ_k}) · σ
```

so the scheme is **unconditionally stable — no Δt ceiling, no sub-stepping**.
For the SLE-coupled / fast-load case the same closed form with the trapezoidal
(linear-in-time σ) variant is used, giving 2nd order in the load without a
stability constraint. (This is *not* the rejected ETD-on-FE-memory route: that
failed because `exp(−M)` per element in the coupled FE basis does not decouple.
Diagonalizing **first** is precisely the decoupling that makes the exponential
exact and cheap.)

## 4. Lateral viscosity (the "LV")

The model is built to run with laterally-varying viscosity. Lateral η breaks the
clean per-degree diagonalization (it couples harmonics), so it is handled the way
FastIsostasy handles its single-mode case, generalized to `K` modes:

1. **Modal basis from a reference profile.** Solve once using the
   **laterally-averaged** radial viscosity profile → `{τ_k(l), residues, b_k}`.
   Lateral variation does **not** trigger a re-solve.
2. **Lateral variation as pointwise rate modulation (split-operator).** Each
   mode's relaxation rate is modulated by the local viscosity on the Gauss grid:
   multiply by the local rate factor in **real space**, apply the degree/mode
   kernel in **spectral space**, transform back with `fe_sht`. This reuses the
   existing 3D-path pattern but acts on `K` **scalar** fields, not a 6-component
   memory tensor — much cheaper.
3. **Depth-weighted per-mode modulation.** Each mode `k` is modulated by the
   lateral viscosity in the **depth band that mode samples**, weighted by the
   mode's radial strain-energy profile (from its mode shape). M0 is driven by
   deep-mantle structure, transient modes by their own layers — a fidelity
   single-mode LV-ELVA cannot represent, and a direct benefit of the multi-mode
   form.

**Exactness ladder (be honest about it):**

- *Radial* η → `K = all` (all modes above tolerance) **converges to** the full
  FE model; `RESP_VE` itself stays available for the exact solution (§7).
- *Lateral* η → modal is an **approximation at every `K`** (true lateral coupling
  is richer than rate modulation). The full tensor-SH path (`RESP_VE`) remains
  the ground-truth reference; modal-LV is the fast, tunable approximation to it,
  cheaper than tensor-SH at all `K`.

## 5. Solving for the modes (matrix-free, dependency-free)

We need the per-degree relaxation spectrum `{τ_k, residues, b_k}`, not a linear
solve — so neither `fe_band` nor LIS applies directly, and we **avoid adding
LAPACK**. We also need only the ~10–30 dominant (slowest-decaying,
surface-coupling) modes, not all ~600. That is the setting for **block subspace
iteration** on the existing, validated relaxation propagator:

- The per-degree relaxation step with the load off —
  `dissipative_rhs → radial_operator_solve_vec → advance_memory` (the
  `fe_advance` path) — IS the discrete propagator `P` applied to a memory vector
  `τ`. `fe_modal` wraps this as a black-box `apply_propagator(l, τ) → Pτ`; no new
  physics, no operator internals touched.
- Iterate a small block (size `K + pad`) of memory vectors through `P`,
  re-orthonormalize each pass, until the subspace stabilizes. The slowest modes
  (`|eig(P)| → 1`) — the physically dominant ones — converge first and fastest.
- The **only** explicit eigensolve is the **Rayleigh–Ritz** projection onto that
  subspace: a tiny `(K+pad)×(K+pad)` matrix (≤ ~30), solved in-house (reduced to
  symmetric via the energy inner product → Jacobi). From it: `e^{Δt s_k} → τ_k`,
  the memory mode shapes (→ residues `r^{u,N,V}_k` via one solve of
  `K⁻¹ D v_k`), and the load projections `b_k`.

The state space is `NLAM × (#Maxwell elements)` ≈ 500–600 per degree, but is
**never formed densely**. Cost is iteration-to-convergence (geometric, fast given
the wide separation of relaxation times), parallel over the `lmax+1` degrees,
**once at init**.

**Fallbacks** if subspace convergence is poor on some degrees: in-house
**Jacobi** on the densely-built `A` (still no dependency), with LAPACK `dgeev`
only as a last resort.

**Runtime per step** is then negligible: `K` exact-exponential scalar updates per
(l,m), versus the per-step banded radial solve of `RESP_VE`. Lateral adds `K`
scalar forward/inverse SHTs per step — cheaper than the 6-component tensor-SH
advance. Net: **faster runtime than full VE, modest one-time init.**

## 6. Architecture / integration

A new response kind `RESP_MODAL` in `fe_response.f90`, beside `RESP_VE`.
Everything downstream is untouched because it goes through the abstract
`response` interface (`response_apply`, `response_horizontal`,
`response_begin/commit_step`, the endpoint brackets, `response_set_dt`,
save/restore, restart).

- **Init** (`response_init_modal`): assemble per-degree operators (existing),
  solve for the modes (§5), select by `mode_rank`, store
  `{τ_k(l), r^{u,N,V}_k(l), b_k(l), weight}`. Reference profile = laterally
  averaged η.
- **State**: `ξ_k,lm` — `K` complex scalars per deforming (l,m), degree-grouped
  like the existing `k` layout. Replaces the memory-stress arrays.
- **apply / horizontal**: affine — elastic gain · σ + Σ_k residue · ξ_k.
- **commit_step / endpoint brackets**: advance `ξ_k` by the exact exponential
  (constant-σ) or its trapezoidal (linear-σ) form; SLE co-convergence uses the
  same closed form re-evaluated against the converged σ.
- **set_dt**: only rescales `e^{−Δt/τ_k}` — no operator re-factor (as today).
- **Lateral** (`response_enable_lateral_visc*`): scalar split-operator via
  `fe_sht` with depth-weighting (§4).
- **Unchanged**: SLE (`fe_sle`), geoid (the per-degree modal sum feeds `gn`),
  rotation/TPW (`fe_rotation`), coupling API (`fe_coupling`), I/O/restart
  (`fe_io`) — restart persists the `ξ_k` fields + the modal basis metadata.
- **Namelist** (`fe_params`, `&fe3d`): `n_modes` (integer; `-1`/`all` = all modes
  above tolerance), `mode_rank` (§2), `mode_residue_tol` (drop modes below this
  fraction of the top residue), and the response/scheme selector.

## 7. Coexistence with `RESP_VE`

`RESP_MODAL` is an **additional** response kind selected by namelist, **not** a
replacement. The full Martinec time-domain solve (`RESP_VE`) is **retained
permanently** and remains selectable as a parameter choice — it is both the exact
solution when wanted and the ground-truth reference for lateral-η validation.
`RESP_MODAL` with `K = all` converges to `RESP_VE` within a set tolerance, but
does not remove or supersede it. The original solving machinery stays in the
code.

## 8. Validation

1. **`K = all`, radial η** converges to the existing FE results within a set
   tolerance ε (M3-L70-V01 Love numbers; 1-D deglaciation; the §3c benchmarks).
2. **Mode-count sweep** of `K` (and `mode_rank`) against full tensor-SH lateral
   runs (M3-L70-V01, deglaciation) → an accuracy/cost curve quantifying the
   approximation per `K`. This curve is the deliverable that justifies "hone to
   context".
3. **Limit checks**: elastic (`t→0⁺`) and isostatic (`t→∞`) limits from the modal
   sum match the directly-computed elastic and fluid Love numbers.

## 9. Implementation plan (incremental, small commits)

**Step 1 — radial-only modal response, validated against `RESP_VE`:**

1. `fe_modal.f90` — `apply_propagator(l,·)` from the existing relaxation kernels,
   block subspace iteration, small in-house Rayleigh–Ritz solve →
   `{τ_k, r^{u,N,V}_k, b_k}`. Pure numerics, unit-testable in isolation.
2. `RESP_MODAL` in `fe_response.f90` — modal state `ξ_k`, `response_init_modal`,
   modal `apply`/`horizontal`/`commit`/endpoint, dispatch wiring; `n_modes`,
   `mode_rank` namelist in `fe_params`.
3. Validate `K = all` converges to `RESP_VE` within tolerance on the 1-D
   benchmark.

**Later steps:**

4. Lateral η via depth-weighted split-operator rate modulation.
5. SLE / endpoint co-convergence + restart wiring for the modal scheme.
6. Accuracy/cost sweep (mode count, `mode_rank`) + docs.
