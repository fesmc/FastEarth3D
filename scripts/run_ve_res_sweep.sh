#!/usr/bin/env bash
#
# run_ve_res_sweep.sh — stage/submit a resolution + sub-stepping sweep of the FULL
# viscoelastic (VE) solver on the real 3-D deglaciation (Bagge 2021 lateral
# viscosity, Tarasov deglaciation forcing, LGM->present). Sister script to
# run_modal_vs_ve.sh, but here there is no modal solver: every run is earth_response=ve
# and we vary only the spherical-harmonic resolution (lmax) and the explicit-scheme
# sub-step ceiling (cfl), measuring how the answer and the cost change.
#
# The production target is lmax=128 (LMAX_REF): high enough for the margins, low
# enough to be affordable. The lmax<128 runs show how far resolution can drop (cost
# falls ~ lmax^2-3) before the answer drifts past tolerance; the cfl probes at the
# reference resolution show whether the default cfl=1 sub-stepping is temporally
# converged. ALL errors are measured against the LMAX_REF / cfl=1 run by
# analysis/ve_res_sweep.jl (it regrids each run onto the reference grid first).
#
# Run it from anywhere; it cd's to the repo root so runme finds .runme/.
#
#   ./scripts/run_ve_res_sweep.sh
#
# Override any setting from the environment, e.g.
#   LMAX_LIST="64 128" ./scripts/run_ve_res_sweep.sh            # fewer resolutions
#   CFL_PROBE="0.5 1.5" ./scripts/run_ve_res_sweep.sh           # extra cfl points at the ref
#   FORCING=/work/ice.nc VISC3D=/work/bagge.nc RUNME_FLAGS="-s -r" ./scripts/run_ve_res_sweep.sh
#   RUNME_FLAGS="-rs -q 12h -w 06:00:00" ./scripts/run_ve_res_sweep.sh
#
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root (where .runme/ lives)

# ============================================================================
# MACHINE-SPECIFIC paths — EDIT THESE for the target cluster (absolute paths).
# ============================================================================
CLIMBER_ROOT=/albedo/work/projects/p_forclima/robinson/models/climber-x
ISOSTASY_DATA=/albedo/work/projects/p_forclima/isostasy_data
FORCING=${FORCING:-${CLIMBER_ROOT}/input/geo_ice_tarasov_deglac.nc}         # ice_thickness(lon,lat,time)
VISC3D=${VISC3D:-${ISOSTASY_DATA}/earth_structure/viscosity/bagge2021.nc}   # log10(eta)(lon,lat,r)

# ============================================================================
# Experiment knobs.
# ============================================================================
# Resolution sweep. LMAX_REF is the production resolution and the error reference;
# it MUST appear in LMAX_LIST. Each run's transform grid is derived from its OWN lmax
# (nlat=2·lmax+2, nphi=4·lmax), so the solver internals run at the resolution set and
# cost falls ~ lmax^2-3 with resolution. Every run uses the one canonical reference
# data/reference/rtopo_gauss_l128.nc (i_eq=1), remapped to its lmax grid online
# (cached) — the input topography stays fixed at l128 regardless of run resolution.
# out.nc thus land on per-run grids; analysis/ve_res_sweep.jl regrids each onto the
# LMAX_REF grid before comparing.
LMAX_LIST=${LMAX_LIST:-32 64 96 128}        # spherical-harmonic degrees to sweep
LMAX_REF=${LMAX_REF:-128}                    # production target + error reference

# Sub-step probe at the reference resolution: extra cfl values (the Maxwell-number
# ceiling M=μΔt/η of the explicit fe scheme; n_sub = ceil(span·max(μ/η)/cfl)). The
# baseline cfl=1.0 IS the LMAX_REF reference run; these add finer (cfl<1, an
# accuracy/convergence check) and coarser (cfl>1, fewer sub-steps) points. Set empty
# to skip. The explicit Maxwell update is stable for M≤2, so cfl=2.0 sits at the
# theoretical edge and cfl=2.5 is expected to break down (NaN) — included on purpose
# to bracket the limit. The analysis tolerates blown-up runs (NaN errors).
CFL_PROBE=${CFL_PROBE:-0.5 1.5 2.0 2.5}

# 3-D-split probe at the reference resolution: extra visc3d_tol values [dex]. An
# element advances on the expensive dyadic SHT path only if its lateral log10(eta)
# spread exceeds visc3d_tol; raising it demotes weakly-3-D elements to their lateral-
# mean 1-D rate, cutting ne3d and the memory-advance cost (~90% of the total) at the
# price of approximating mild lateral structure. Default (1e-3 dex) IS the reference;
# these add coarser splits. Set empty to skip.
VTOL_PROBE=${VTOL_PROBE:-0.1 0.3 1.0}

T0=${T0:--26000.0}                           # transient start [yr] (LGM)
T1=${T1:-0.0}                                # transient end   [yr] (present)
DT_COUPLE=${DT_COUPLE:-100.0}                # coupling interval [yr] (forcing cadence)
EQUIL_TIME_MAX=${EQUIL_TIME_MAX:-10000.0}    # LGM-memory spin-up cap [yr]
OMP=${OMP:-8}                                # OpenMP threads per run
EXP=${EXP:-runs/ve_res_sweep}                # experiment root (under gitignored runs/)

# How runme launches each run (see run_modal_vs_ve.sh):
#   "-s -r" prepare SLURM scripts AND submit (HPC, default); "-s" stage only;
#   "-r" run locally in the background; "" stage the run dirs only.
RUNME_FLAGS=${RUNME_FLAGS--s -r}             # note: `-` not `:-`, so RUNME_FLAGS="" = stage-only

# 3-D viscosity + cheap 1-D LGM spin-up (the deglac3d setup from run_modal_vs_ve.sh).
VISC3D_ON=(fe3d.l_visc_3d=true fe3d.pre_spinup_1d=true fe3d.visc_3d_file="$VISC3D")

# Parameters common to every run (machine paths + the shared deglaciation setup).
# Resolution (lmax + reference) and the scheme/cfl are added per run below; the grid
# is derived from each run's own lmax (not pinned).
COMMON=(
  fe3d.file_forcing="$FORCING"
  fe3d.name_ice=ice_thickness
  fe3d.i_eq=1
  fe3d.earth_response=ve
  fe3d.scheme=fe
  fe3d.dt_couple="$DT_COUPLE"
  fe3d.equil_time_max="$EQUIL_TIME_MAX"
  fe3d.time_init="$T0"
  fe3d.time_end="$T1"
  fe3d.rotation=true            # real-Earth run: rotational feedback on
  fe3d.file_out=out.nc
)

# the single canonical reference; remapped to each run's resolution online (cached).
REF=${REF:-data/reference/rtopo_gauss_l128.nc}
[ -f "$REF" ] || { echo "ERROR: reference file not found: $REF (make fastearth_mkref)" >&2; exit 1; }

# launch <outdir> <lmax> <extra -p args...> — one VE run at the given resolution.
launch() {
  local out=$1 lmax=$2; shift 2
  local ref="$REF"
  echo ">>> $out   (lmax=$lmax $*)"
  runme -o "$out" -e main --omp "$OMP" $RUNME_FLAGS \
        -p "${COMMON[@]}" \
        fe3d.lmax="$lmax" fe3d.z_bed_ref_file="$ref" fe3d.h_ice_ref_file="$ref" \
        "${VISC3D_ON[@]}" "$@"
}

# LMAX_REF must be one of the swept resolutions (it is the error reference).
case " $LMAX_LIST " in *" $LMAX_REF "*) ;; *) echo "ERROR: LMAX_REF=$LMAX_REF not in LMAX_LIST='$LMAX_LIST'" >&2; exit 1;; esac

echo "lmax sweep: $LMAX_LIST   ref=$LMAX_REF   cfl-probe@ref: ${CFL_PROBE:-(none)}   vtol-probe@ref: ${VTOL_PROBE:-(none)}"
echo "window=[$T0,$T1]  dt_couple=$DT_COUPLE  equil_time_max=$EQUIL_TIME_MAX  omp=$OMP"
echo "exp root: $EXP    runme flags: '$RUNME_FLAGS'"

# ---------------------------------------------------------------------------
# Resolution sweep (explicit fe scheme, cfl=1). lmax.<L>; lmax.<REF> is the reference.
# Each run's grid is derived from its own lmax (see COMMON note above).
# ---------------------------------------------------------------------------
for L in $LMAX_LIST; do
  launch "$EXP/lmax.$L" "$L" fe3d.cfl=1.0
done

# ---------------------------------------------------------------------------
# Sub-step probe at the reference resolution: extra cfl values (cfl=1 already done).
# ---------------------------------------------------------------------------
for C in ${CFL_PROBE:-}; do
  launch "$EXP/lmax.$LMAX_REF.cfl$C" "$LMAX_REF" fe3d.cfl="$C"
done

# ---------------------------------------------------------------------------
# 3-D-split probe at the reference resolution: extra visc3d_tol values (default done).
# ---------------------------------------------------------------------------
for V in ${VTOL_PROBE:-}; do
  launch "$EXP/lmax.$LMAX_REF.vtol$V" "$LMAX_REF" fe3d.cfl=1.0 fe3d.visc3d_tol="$V"
done

echo "done."
