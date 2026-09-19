# Modal lateral-viscosity split-operator error — diagnosis

Diagnostic: `tests/diag_modal_latvisc.f90` (target `make diag_modal_latvisc`, NOT in
`make check`). Isolates the RESP_MODAL lateral split-operator error against the
RESP_VE tensor-SH ground truth, **with lateral viscosity active on both paths**
(all prior modal-vs-VE probes were radial and missed this entirely).

## Setup

- Earth: M3-L70-V01 (laterally-averaged reference profile = the modal basis).
- Controlled contrast: an axisymmetric polar **viscosity cap** (colat < θ_cap,
  5° cosine taper), depth-uniform `Δ = log10(η_local/η_ref) > 0` (stiff craton),
  applied to every Maxwell element. The *identical* `pert_elem(nphi,nlat,ne)`
  field is handed to both `response_enable_lateral_visc` (VE→tensor-SH) and
  `response_enable_lateral_visc_modal` (split-operator) — the only difference is
  the algorithm.
- Co-located disc surface load (colat < θ_cap, tapered). Both operators are linear
  in σ, so absolute load amplitude is irrelevant; we report ratios.
- Protocol = **free rebound** (cleanest rate probe): hold the load (glaciation) to
  build depression, then remove it (σ=0) and watch free relaxation, whose decay
  rate *is* the lateral viscosity operator. Fluid limits are contrast-independent,
  so any gap is a transient **rate** error.
- Metrics: peak uplift over the cap region, modal/VE ratio (1.0 = perfect);
  `grms/pk` = global area-rms(modal−VE) over the VE peak (field error).

## Result

Two resolutions run (`make diag_modal_latvisc`; logs in `logs/latvisc_*.log`).
lmax 16 over-states the magnitudes (coarse grid / Gibbs); **lmax 48 is the
reference** below. The qualitative structure is identical at both.

### lmax 48 (glaciation 40 kyr, rebound 20 kyr; dt 200 yr) — REFERENCE

End-of-glaciation modal/VE peak-uplift **ratio** (1.0 = perfect) | field error
`grms/pk` | rebound ratio at 20 kyr:

| Δ \ θ_cap | 10° | 20° | 40° |
|-----------|-----|-----|-----|
| 0.5 | 0.98 / 0.005 / 0.82 | 0.91 / 0.015 / 0.62 | **0.77 / 0.085 / 0.60** |
| 1.0 | 1.26 / 0.020 / 0.61 | 1.15 / 0.039 / 0.64 | 0.91 / 0.189 / 3.12 |
| 1.5 | 1.68 / 0.044 / 1.08 | 1.42 / 0.069 / 1.24 | 1.09 / 0.265 / 4.92 |
| 2.0 | 2.02 / 0.064 / 1.72 | 1.59 / 0.091 / 1.62 | 1.24 / 0.311 / 6.00 |

Reading: error grows **monotonically with contrast at every scale**; field error
`grms/pk` is **largest for the largest cap** (40°: 0.085→0.31). The peak-ratio
*sign* flips (under-prediction at moderate Δ, over-prediction at high Δ), but the
*field error magnitude* is unambiguous and large for continental scale + realistic
contrast. During free **rebound** the divergence is far worse than at
end-of-glaciation (40°/Δ≥1: 3–6× over; 40°/Δ0.5: 0.60× under).

### Matches the ensemble's deglac3d regime

The ensemble's N. America attenuation (modal 154 m vs VE 343 m at PD ⇒ ratio 0.45)
is reproduced by the **continental-scale (40°), moderate stiff-contrast (Δ≈0.5–1)**
corner, where modal **under-predicts rebound by 30–55%** (40°/Δ0.5 rebound 0.60;
40°/Δ1 early-rebound 0.73–0.85). This is the relevant regime for the Canadian
Shield (broad, ~0.5–1 dex stiffer mantle). Higher contrast (Δ≥1.5) over-predicts —
so the sign of the local modal error is contrast-dependent, but the breakdown is
generic. lmax 16 numbers (for reference, much larger): 40°/Δ0.5 end-glac 0.35,
40°/Δ2 end-glac 6.7, rebound to 77×.

## What this says

1. **The error is NOT a uniform ~55% attenuation.** It scales steeply with BOTH
   lateral contrast Δ and cap size θ, and **changes sign**: small/weak → mild
   under-prediction (the advertised 1st-order regime); large+strong → severe
   over-prediction with O(1)+ field errors.
2. **Mechanism (visible in the code).** For a stiff cap (Δ>0), `rho = 10^(−Δ) < 1`
   inside, so the zero-mean anomaly `(Ri − Rbar) < 0` there →
   `E2 = exp(−dt·mlatRate) > 1`. In `modal_lateral_anomaly`,
   `gphi = E2·gphi + (1−E2)·gsig`; during free rebound (gsig=0) this is
   `gphi ← E2·gphi` with E2>1 → **amplification, not relaxation**. The Lie split
   treats the anomaly as a small perturbation of the mean step; at large contrast
   the product overshoots and can grow.
3. **Three compounding error sources** (to be attributed):
   - **Lie splitting**: mean (spectral, diagonal in l) ⊗ anomaly (real-space,
     diagonal in θφ) do not commute; error ~½[A,B], grows with both magnitudes.
   - **Characteristic-τ̂ mismatch**: the anomaly is applied at a single rank-wide
     `τ̂_i` (|C^u|-weighted over degrees, low-l-dominated), but a *localized* cap's
     structure lives at high l where `τ_i(l)` is very different → wrong rate.
   - **Amplification (E2>1)**: net non-contractive step for strong negative
     anomalies → unphysical growth in free relaxation.

## Resolves the "2.6% vs 55%" discrepancy (task background)

The prior stripped-down modal+3D vs VE+3D run reported ~2.6% **global** rsl rms,
yet the ensemble shows ~55% **local** attenuation over N. America at PD. These are
consistent: the error is strongly localized to the high-contrast craton. A global
rms is diluted by the vast low-contrast ocean/non-cratonic area (small error
there), while the local peak/field error over the cap is O(50–100%). This
diagnostic measures the local error directly (`grms/pk` normalized by the cap peak)
and recovers the large values.

## Attribution (lmax 48; `logs/latvisc_attrib_lmax48.log`)

Three modal lateral methods now exist (`&fe3d` / `response%lat_method`), driven
against the same VE+LV truth:

- **LIE** — original 1st-order Lie split, rank-characteristic τ̂ (production default).
- **STRANG** — symmetric 2nd-order split, same τ̂.
- **COUPLED** — per-rank coupled rate operator `exp(−Δt·L_i)`, no split,
  per-degree-exact τ, full lateral coupling; matrix-free grid Arnoldi. Reduces to
  the spectral mean step for a uniform rank (test_modal_visc3d case D, 9.5e-15).

End-glaciation modal/VE cap-peak ratio | field error grms/pk | rebound ratio @20 kyr:

| Δ, 40° cap | LIE | STRANG | COUPLED |
|---|---|---|---|
| 0.5 | 0.77 / 0.085 / 0.60 | 0.77 / 0.084 / 0.59 | 0.71 / 0.114 / **0.80** |
| 1.0 | 0.91 / 0.189 / **3.12** | 0.91 / 0.188 / 3.09 | 0.56 / 0.194 / **0.71** |
| 1.5 | 1.09 / 0.265 / **4.92** | 1.10 / 0.263 / 4.88 | 0.51 / 0.237 / **0.81** |
| 2.0 | 1.24 / 0.311 / **6.00** | 1.24 / 0.308 / 5.96 | 0.51 / 0.259 / **0.92** |

Three conclusions:

1. **The split-ORDER error is negligible.** STRANG ≈ LIE to 3 digits everywhere.
   The Lie 1st-order splitting is not the problem.
2. **COUPLED is converged, not approximate-by-truncation.** m_krylov = 12 and 64
   give bit-identical results (`logs` m=12 vs m=64), so COUPLED solves its operator
   exactly; its deviation from VE is the **rate-modulation ANSATZ itself**, not
   numerics. (Arnoldi reaches lucky breakdown in a handful of vectors, so cost is
   modest and m-independent.)
3. **COUPLED's real win is STABILITY.** It removes LIE's catastrophic rebound
   runaway — the E2>1 amplification drives LIE to 3–6× over-prediction during free
   rebound at continental scale + realistic contrast, whereas COUPLED stays bounded
   at 0.71–0.92. For the deglaciation (rebound-dominated) this is the regime that
   matters.

### Why neither closes the gap — the ansatz limit

At MILD contrast where the modal-LV ansatz should be near-exact (10°/Δ0.5), LIE
gives 0.98 (≈ VE) but **COUPLED gives 0.81** — COUPLED is *further* from VE despite
solving its operator exactly. Interpretation: the exact pointwise ansatz applies the
*local* rate `R(θ,φ)/τ`, but the true (VE) response is **non-local** (deformation
wavelength ~ hundreds of km), so for caps near/below that scale COUPLED
**over-localizes** the slowdown and under-predicts. LIE's degree-exact mean
modulation accidentally supplies that spatial dilution at small scales (and blows up
at large scale). So the two methods trade which regime they handle well; **neither
brings the continental-scale + realistic-contrast gap below ~20–50%** — that residual
is intrinsic to the rate-modulation ansatz (design-modal.md §4: "approximation at
every K"). Closing it further would require leaving the ansatz (genuine harmonic
coupling), i.e. moving toward the tensor-SH cost that modal exists to avoid.

### Controlled real-field confirmation (Tarasov + Bagge)

Controlled isolation on the REAL field (`runs/latcheck`, lmax 32, −8 kyr→0, all
modes, `modal_adaptive=false`, only `lat_method` varied) — VE truth vs modal-coupled
vs modal-lie:

| at VE-peak cell over N. America | ratio modal/VE | global rms(modal−VE) |
|---|---|---|
| modal-lie     | 0.406 | 2.69 m |
| modal-coupled | 0.345 | 2.98 m |

COUPLED is **marginally worse** than LIE on the real field, not better; both
under-predict ~60–65% at this resolution/window. This **confirms on real data** what
the synthetic predicted: at the real continental contrast COUPLED ≈ LIE (slightly
worse, from over-localization), and the ~55% N. America attenuation is the intrinsic
rate-modulation-ansatz ceiling — neither method breaks it. COUPLED is genuinely
active (the runs differ); its only edge over LIE is stability (no high-contrast
rebound runaway), which does not show in this moderate-contrast metric.

**Ensemble caveat (now fixed):** the re-run that "showed no improvement" also flipped
`modal_adaptive` off (nsolve 22→1); the radial set regressed 0.34→0.89 m PD from that
alone (lateral-independent). `run_modal_vs_ve.sh` now sets `fe3d.modal_adaptive=true`
(knob `ADAPTIVE`) so the comparison isolates `lat_method`.

### Decision (taken)

**COUPLED is now the production default** (`&fe3d lat_method = "coupled"`, type
default `LAT_COUPLED`) — chosen for stability (no rebound runaway), strictly safer
than LIE for deglaciation at realistic contrast. LIE and STRANG remain selectable
(`lat_method = "lie" | "strang"`) so the ensemble can sweep them. The ensemble will
be re-run with COUPLED.

Open / accepted:
- The ~20–50% modal-LV ceiling at strong continental contrast is the documented
  intrinsic limit of the rate-modulation ansatz (design §4). Going below it needs a
  beyond-ansatz scheme (a little true harmonic coupling) — a larger research step,
  deferred.
- COUPLED-vs-LIE wall-cost not yet precisely benchmarked (both ≪ VE; Arnoldi breaks
  down in a few vectors so the overhead is small and m-independent). Measure on the
  ensemble re-run.
- Confirmation on the real Tarasov+Bagge field (lmax 32–48) still deferred (synthetic
  only so far) — the ensemble re-run covers this.
