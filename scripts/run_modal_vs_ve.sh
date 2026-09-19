#!/usr/bin/env bash
#
# run_modal_vs_ve.sh — stage/submit the experiment set that quantifies how well
# the reduced MODAL response approximates the full viscoelastic (VE) solver, in
# both accuracy and cost, using runme (https://github.com/.../runme).
#
# Two sets, both forced by the Tarasov deglaciation and measured against the VE
# run in the *same* set (the ground truth):
#
#   1. radial   — 1-D (radially-symmetric) viscosity. The clean limit where modal
#                 with n_modes=all converges to VE exactly; isolates the accuracy
#                 of the mode-count / ranking dial (no lateral approximation).
#   2. deglac3d — full run with laterally-varying (3-D) viscosity. The LGM spin-up
#                 is done with the cheap 1-D solver in every case (pre_spinup_1d=true),
#                 then the transient runs on the 3-D path. This is where modal is a
#                 genuine approximation to VE (design §4).
#
# Each set sweeps the modal dial: earth_response=modal x n_modes x mode_rank, plus
# n_modes=all, plus the single VE reference. Compare each modal out.nc (rsl/bsl)
# against the set's ve/out.nc, and the wall time from out.out ([PROFILE] lines) or
# the SLURM accounting.
#
# Run it from anywhere; it cd's to the repo root so runme finds .runme/.
#
#   ./scripts/run_modal_vs_ve.sh
#
# Override any setting from the environment, e.g.
#   LMAX=128 ./scripts/run_modal_vs_ve.sh                      # whole script at lmax 128
#   LMAX_DEGLAC=128 ./scripts/run_modal_vs_ve.sh               # radial at 64, deglaciation at 128
#   FORCING=/work/ice.nc VISC3D=/work/bagge.nc RUNME_FLAGS="-s -r" ./scripts/run_modal_vs_ve.sh
#   RUNME_FLAGS="-rs -q 12h -w 06:00:00" ./scripts/run_modal_vs_ve.sh 
#
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root (where .runme/ lives)

# ============================================================================
# MACHINE-SPECIFIC paths — EDIT THESE for the target cluster (absolute paths).
# ============================================================================
CLIMBER_ROOT=/albedo/work/projects/p_forclima/robinson/models/climber-x
ISOSTASY_DATA=/albedo/work/projects/p_forclima/isostasy_data
ICE_DATA=/albedo/work/projects/p_forclima/ice_data
FORCING=${FORCING:-${CLIMBER_ROOT}/input/geo_ice_tarasov_deglac.nc}   # ice_thickness(lon,lat,time)
VISC3D=${VISC3D:-${ISOSTASY_DATA}/earth_structure/viscosity/bagge2021.nc}  # log10(eta)(lon,lat,r)

# ============================================================================
# Experiment knobs.
# ============================================================================
# Resolution. LMAX is the one knob for the WHOLE script; LMAX_RADIAL / LMAX_DEGLAC
# default to it but can be set independently (e.g. cheap radial benchmark at 64,
# deglaciation at 128). All runs use the canonical reference rtopo_gauss_l128.nc.
LMAX=${LMAX:-64}                             # spherical-harmonic degree (whole script)
LMAX_RADIAL=${LMAX_RADIAL:-$LMAX}            # set 1 (idealized / radial)
LMAX_DEGLAC=${LMAX_DEGLAC:-$LMAX}            # set 2 (full deglaciation)
T0=${T0:--26000.0}                           # transient start [yr] (LGM)
T1=${T1:-0.0}                                # transient end   [yr] (present)
DT_COUPLE=${DT_COUPLE:-100.0}                # coupling interval [yr] (forcing cadence)
EQUIL_TIME_MAX=${EQUIL_TIME_MAX:-100000.0}    # LGM-memory spin-up cap [yr]
OMP=${OMP:-8}                                # OpenMP threads per run
EXP=${EXP:-runs/modal_vs_ve}                 # experiment root (under gitignored runs/)

# Modal dial swept in each set (comma lists => runme ensemble dimensions).
N_MODES=${N_MODES:-1,2,4,8}                  # modes kept per degree (truncated)
RANKS=${RANKS:-isostatic,rate,residue}       # mode_rank metric

# Adaptive step-doubling for the MODAL response. modal_adaptive defaults OFF in the
# model (1 exact-exp step/couple — fast, but a temporal-truncation error that A3
# removes), so the ensemble must turn it ON to sub-step; otherwise RTOL is inert and
# modal runs 1-step (radial accuracy degrades, deglac3d unaffected). RTOL is the
# step-doubling tolerance: ~1e-4 best accuracy (~9 SLE/couple), 1e-3 a balance
# (~3.6 SLE/couple). VE uses scheme=fe and ignores both. See diag_modal_sle.
ADAPTIVE=${ADAPTIVE:-true}        # modal: sub-step each couple to RTOL (A3). false = 1-step
RTOL=${RTOL:-1e-4}

# How runme launches each (group) of runs:
#   "-s -r"  prepare SLURM scripts AND submit   (HPC, default)
#   "-s"     prepare SLURM scripts, do NOT submit (inspect, submit by hand)
#   "-r"     run locally in the background
#   ""       stage the run dirs only (no run)
RUNME_FLAGS=${RUNME_FLAGS--s -r}     # note: `-` not `:-`, so RUNME_FLAGS="" means stage-only

# Turning on 3-D viscosity + 1-D spin-up. runme writes booleans quoted ('true'),
# but the model's nml reader parses 'true'/'false' as logicals, so -p is fine.
VISC3D_ON=(fe3d.l_visc_3d=true fe3d.pre_spinup_1d=true fe3d.visc_3d_file="$VISC3D")

# Per-resolution params (lmax, i_eq=1). Every run uses the one canonical present-day
# reference, remapped to its lmax online (cached). The ref is resolved via the rundir
# 'data' symlink; require it to exist up front.
REF=${REF:-data/reference/rtopo_gauss_l128.nc}
[ -f "$REF" ] || { echo "ERROR: reference file not found: $REF (make fastearth_mkref)" >&2; exit 1; }
RES_RADIAL=(fe3d.lmax="$LMAX_RADIAL" fe3d.z_bed_ref_file="$REF" fe3d.h_ice_ref_file="$REF")
RES_DEGLAC=(fe3d.lmax="$LMAX_DEGLAC" fe3d.z_bed_ref_file="$REF" fe3d.h_ice_ref_file="$REF")

# Parameters common to every run (machine paths + the shared deglaciation setup).
# Resolution (lmax + reference) is added per set from RES_RADIAL / RES_DEGLAC.
COMMON=(
  fe3d.file_forcing="$FORCING"
  fe3d.name_ice=ice_thickness
  fe3d.i_eq=1
  fe3d.dt_couple="$DT_COUPLE"
  fe3d.equil_time_max="$EQUIL_TIME_MAX"
  fe3d.time_init="$T0"
  fe3d.time_end="$T1"
  fe3d.rotation=true            # real-Earth runs: rotational feedback on (both solvers)
  fe3d.modal_adaptive="$ADAPTIVE" # modal: enable A3 sub-stepping (VE ignores it)
  fe3d.rtol="$RTOL"             # modal adaptive sub-step tolerance (VE ignores it)
  fe3d.file_out=out.nc
)

# launch <outdir> <ensemble:0|1> <extra -p args...>
launch() {
  local out=$1 ens=$2; shift 2
  local aflag=()
  [ "$ens" = "1" ] && aflag=(-a)                     # name ensemble member dirs from their params
  echo ">>> $out   ($*)"
  runme -o "$out" -e main --omp "$OMP" $RUNME_FLAGS \
        ${aflag[@]+"${aflag[@]}"} \
        -p "${COMMON[@]}" "$@"
}

echo "lmax: radial=$LMAX_RADIAL deglac3d=$LMAX_DEGLAC  window=[$T0,$T1]  dt_couple=$DT_COUPLE  equil_time_max=$EQUIL_TIME_MAX  omp=$OMP"
echo "exp root: $EXP    runme flags: '$RUNME_FLAGS'"

# ---------------------------------------------------------------------------
# Set 1 — idealized: 1-D (radial) viscosity. modal(n_modes=all) -> VE exactly.
# ---------------------------------------------------------------------------
launch "$EXP/radial/ve"        0 "${RES_RADIAL[@]}" fe3d.earth_response=ve    fe3d.scheme=fe          fe3d.l_visc_3d=false
launch "$EXP/radial/modal"     1 "${RES_RADIAL[@]}" fe3d.earth_response=modal fe3d.n_modes="$N_MODES" fe3d.mode_rank="$RANKS" fe3d.l_visc_3d=false
launch "$EXP/radial/modal_all" 0 "${RES_RADIAL[@]}" fe3d.earth_response=modal fe3d.n_modes=-1         fe3d.l_visc_3d=false

# ---------------------------------------------------------------------------
# Set 2 — full deglaciation: 3-D viscosity, 1-D spin-up.
# ---------------------------------------------------------------------------
launch "$EXP/deglac3d/ve"        0 "${RES_DEGLAC[@]}" fe3d.earth_response=ve    fe3d.scheme=fe          "${VISC3D_ON[@]}"
launch "$EXP/deglac3d/modal"     1 "${RES_DEGLAC[@]}" fe3d.earth_response=modal fe3d.n_modes="$N_MODES" fe3d.mode_rank="$RANKS" "${VISC3D_ON[@]}"
launch "$EXP/deglac3d/modal_all" 0 "${RES_DEGLAC[@]}" fe3d.earth_response=modal fe3d.n_modes=-1         "${VISC3D_ON[@]}"

echo "done."
