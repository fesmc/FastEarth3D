# The optional VILMA backend (`solver = "vilma"`)

FastEarth3D can run with **VILMA** (Martinec/Klemann — the CLIMBER-X `i_geo=2`
solid-earth backend) in place of its own solver, behind the *same* driver,
namelist, forcing, remap and output. That makes an F-vs-V comparison a one-line
namelist change instead of two different experiments: identical ice forcing,
identical reference state, identical output grid and file layout.

This is a **core-developer facility**, not a user feature.

---

## 1. VILMA is NOT a dependency

The default build (`make fastearth`) does not reference a single VILMA symbol,
needs no `vilma/include`, and does not link `vega_pism.a`. VILMA is a
hand-installed, precompiled library that is absent on most machines, and it must
stay that way.

Everything VILMA-specific lives inside `#ifdef VILMA` in `src/fe_vilma.f90`,
which the `vilma=1` make switch turns on (`config/common.mk`, mirroring
CLIMBER-X's own `vilma=` / `fastearth=` toggles). With `vilma=0` — the default —
that file compiles to a pure-Fortran stub, and asking for `solver="vilma"` at
runtime aborts at init with an actionable message:

```
 ======================================================================
  FastEarth3D: solver="vilma" requested, but this binary has no VILMA.
 ======================================================================
  ...
  To use it, rebuild with the backend switched on:

      make clean
      make fastearth vilma=1 VILMAROOT=/path/to/vilma
  ...
```

Not a link error, not a crash, not a silent no-op.

---

## 2. Installing VILMA

VILMA is distributed as a precompiled library, not as source you build here:

```
$VILMAROOT/
  include/*.mod      # Intel .mod files (mod_firstlevel, mod_sle, mod_io, ...)
  lib/vega_pism.a    # the library archive
```

It is cloned/installed by hand (at DKRZ: `git clone git@github.com:cxesmc/vilma.git`,
which ships the prebuilt `include/` and `lib/`). It also needs its own static
input files — `densi.inp`, `tint.inp`, `visko.inp` (or a 3-D viscosity NetCDF)
and a `vilma_grid.nc` describing its own grid. The CLIMBER-X checkout carries a
working set under `input/vilma/`, which is what `vilma_input_dir` should point at.

Two hard constraints on the install:

* **Compiler.** The `.mod` files are Intel-specific and version-stamped. Only a
  compatible `ifort`/`ifx` can read them; gfortran cannot.
* **netCDF-Fortran ABI.** See §6 — this is the one that will actually bite you.

VILMA needs no LAPACK, no BLAS and no MPI (it carries its own banded solvers and
FFT), but it *is* OpenMP-parallel, so build with `openmp=1` (the default).

---

## 3. Building and selecting it

```sh
make clean
make fastearth vilma=1 VILMAROOT=/work/ba1442/robinson/models/vilma
```

Then, in the run config's `&fe3d` group:

```fortran
&fe3d
    solver = "vilma"

    lmax   = 32                       ! the COUPLING + OUTPUT grid (see §4)

    vilma_jmax         = 170          ! VILMA's own spectral degree
    vilma_input_dir    = "input/vilma"
    vilma_out_dir      = "vilma"      ! VILMA scratch/output, created if absent
    vilma_grid_file    = "input/vilma/vilma_grid.nc"
    vilma_visc_1d_file = "visko.inp"
    vilma_visc_3d_file = "visc3d_Bagge2021.nc"   ! used only when l_visc_3d = .true.
/
```

`solver` defaults to `"fe3d"`, so **existing namelists are unaffected**.

When `solver="vilma"`, the `&fe3d` settings that describe the *FastEarth3D*
solver — `earth`/layer arrays, `scheme`, `earth_response`, the `sle_*` knobs, the
adaptive-Δt knobs, `rotation`, `equil_*`/`pre_spinup_1d` — are **ignored**. VILMA
has its own earth structure, its own sea-level equation and its own time
stepping. Only `lmax`/`nlat`/`nphi` (the coupling and output grid) and
`l_visc_3d` (which sets VILMA's `vg%l_mod`) are shared.

---

## 4. Which grid is which

Three grids, and the boundaries matter:

| grid | size (example) | who owns it |
|---|---|---|
| host lon-lat | GLAC-1D's native grid | the forcing file; `fe_drive`/`fe_remap` map it away |
| **model Gauss** | `lmax=32` → 128 × 66 | FastEarth3D's Gauss-Legendre grid (`fe_sht`) |
| VILMA grid | `jmax=170` → 512 × 256 | VILMA's own Gauss grid |

* The driver remaps the forcing onto the **model Gauss grid**, exactly as for the
  native solver (`fe_remap`, conservative).
* `fe_vilma` owns the second leg and nothing else:
  * model Gauss → VILMA: **conservative** (`coords` `"con"`), for ice thickness
    and the reference fields;
  * VILMA → model Gauss: **bilinear**, for relative sea level.
  Weights are cached under `maps/`, keyed by grid dimensions.
* `se%gg%*`, `se%rsl`, `se%z_bed` and the whole output file are on the **model
  Gauss grid** in both cases, which is what makes the comparison direct.
* VILMA's `rsl` is `real(dp)(nlat, nlon)` in metres, latitude-first and
  north-first, positive = sea-level rise. `fe_vilma` transposes it before
  remapping, and checks its shape against `vilma_grid_file` every step.

`vilma_jmax` and `vilma_grid_file` must describe the same grid: VILMA derives its
working grid from `jmax` alone (`nlon` = smallest power of two > `3*jmax`,
`nlat = nlon/2`; so `jmax=170` → 512 × 256 = N128). `fe_vilma_init` checks this
and refuses a mismatch with a message naming both settings.

Two deliberate differences from CLIMBER-X's wrapper, both in the direction of a
fairer comparison:

* the ice delivered to VILMA is **globally rescaled** so its mass equals the mass
  on the model Gauss grid as SHTns integrates it (CLIMBER-X does not rescale);
* `loadh.inp` is written with FastEarth3D's **own** `rho_ice`/`rho_water`
  (931/1000) rather than CLIMBER-X's hard-coded 910/1020, so both backends turn
  the same ice thickness into the same load.

---

## 5. What is populated, and what is not

`solid_earth` and `fe_io` expose more than VILMA produces. Nothing here is faked.

### Populated (derived honestly from VILMA output)

| field | how |
|---|---|
| `rsl` | VILMA's `mod_sle::rsl`, remapped to the model Gauss grid |
| `z_bed` | `z_bed_eq - rsl`, the same reconstruction the native solver uses |
| `h_ice` | the forcing, as remapped |
| `C_ocean` | `fe_sle`'s `ocean_function` (flotation) applied to VILMA's updated bed and the current ice. This is **FastEarth3D's diagnostic of VILMA's state**, not VILMA's internal ocean function — the library keeps `ocfunc` private and does not expose it on this grid. The two can differ by a cell where the flotation rules disagree. |
| `bsl` | `update_bsl`'s barystatic integral over that `C_ocean` and the ice anomaly — again a FastEarth3D diagnostic of VILMA's state, computed by the *same* formula for both backends. VILMA's own water-layer scalar is written to `<vilma_out_dir>/vega_oce.dat`, column 4 (columns 5–6 are its degree-0 load and ocean-area fraction) if you want VILMA's internal number instead. |
| `t_solver` | wall-clock inside `time_evolution` |
| `t_remap` | wall-clock in the Gauss ↔ VILMA remap |

### Not populated (documented fill value, writer/printer skips it)

| field | why | value |
|---|---|---|
| `worst_mass_resid` | VILMA reports no sea-level-equation mass residual. | `FE_UNSET = -9999` (`fe_coupling`). The driver omits it from the per-step line instead of printing it; it is not written to NetCDF by either backend. |
| `resp%t_drift`, `resp%t_mem`, `sle%t_*`, `stepper%t_guard`, `stepper%n_*` | those phases do not exist — there is no drift solve, no Maxwell memory advance, no FastEarth3D SLE iteration and no adaptive sub-stepping. | left at 0; `fe_drive` **skips** the two native breakdown blocks and the sub-step line entirely and prints a VILMA-specific breakdown (`time_evolution` / remap / other) instead. |
| `dt_try` and the `tau_*` / `phi_*` memory fields | no FastEarth3D memory state exists. | the response is initialised NULL, so `fe_io` writes only the common diagnostic set and no memory dimensions. |

### Refused rather than approximated

* **`solid_earth_spinup`** (`equil_time_max > 0` or `pre_spinup_1d`) aborts with
  an explanatory message. It relaxes *this model's* viscous memory against a held
  reference; VILMA advances its memory only along its own time axis, and faking
  it would silently consume VILMA's clock. Set `equil_time_max = 0` and
  `pre_spinup_1d = .false.` and drive the full window instead, which is how
  CLIMBER-X runs VILMA.
* **`restart_in_file`** is refused: a FastEarth3D restart carries the native
  solver's memory, which has no VILMA counterpart. VILMA restarts through its own
  `r_restart`/`w_restart` files, which this driver does not wire up — so no
  FastEarth3D restart is written at the end of a VILMA run either.

---

## 6. Known gaps and gotchas

### 6a. netCDF-Fortran ABI — read this before anything else

**One binary, one netCDF-Fortran ABI.** `nf90_def_var` gained trailing optional
arguments (`quantize_mode`, `nsd`) in netCDF-Fortran **4.6.0**, which moved the
hidden character-length argument. Compile against one side and link the other and
the name length is read from the wrong slot:

* compiled ≥4.6, linked 4.5.x → `NetCDF: Name contains illegal characters`
* compiled 4.5.x, linked ≥4.6 → segfault in `nf_def_var_quantize`

The precompiled VILMA at DKRZ was built against **netCDF-Fortran ≥ 4.6**, while
FastEarth3D's `configme`-detected default on levante is **4.5.3**. Linking them
as-is builds and runs all the way through `setup` and most of the first
`time_evolution`, then dies inside `sle_write_rsl`:

```
 GCM :     0.000000E+00     0.000000E+00     0.000000E+00
  => ERROR: NetCDF: Name contains illegal characters
  STOP
```

Minimal reproducer (`nf90_def_dim` is fine, `nf90_def_var` is not):

```sh
ifx -I/sw/spack-levante/netcdf-fortran-4.6.2-5t6lbs/include dv.f90 -o dv.x \
    -L/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib -lnetcdff ... && ./dv.x
 create No error
 def_dimNo error
 def_varNetCDF: Name contains illegal characters
```

So **everything in the binary — fesm-utils/ncio included — must be built against
a netCDF-Fortran that matches VILMA's.** On levante that means an Intel-built
4.6.x, e.g. `netcdf-fortran-4.6.2-5t6lbs` with `netcdf-c-4.9.2-x7g75q`.

Because `fesm-utils` is a *shared* checkout, do not re-point it in place. Build a
private copy and pass it on the make line — `FESMUTILSROOT`, `INC_NC` and
`LIB_NC` are all overridable:

```sh
# 1. a private fesm-utils against the matching netCDF
cp -r <fesm-utils>/{src,config,Makefile,config.py} $FU/
ln -s <fesm-utils>/{fftw,SHTns,lis} $FU/          # no netCDF in these
#    edit $FU/Makefile: NC_FROOT / NC_CROOT / INC_NC / LIB_NC -> the 4.6.x pair
( cd $FU && make fesmutils-static openmp=1 )      # serial: its -j deps are incomplete

# 2. FastEarth3D against that, plus VILMA
A=/sw/spack-levante/netcdf-fortran-4.6.2-5t6lbs
C=/sw/spack-levante/netcdf-c-4.9.2-x7g75q
make clean
make fastearth vilma=1 VILMAROOT=/work/ba1442/robinson/models/vilma \
     FESMUTILSROOT=$FU \
     INC_NC="-I$A/include -I$C/include" \
     LIB_NC="-L$A/lib -lnetcdff -L$C/lib -lnetcdf -Wl,-rpath,$A/lib -Wl,-rpath,$C/lib"
```

This is an environment constraint, not a code one, so it is deliberately *not*
baked into `config/common.mk`: the default build keeps the machine's configured
netCDF untouched. The clean long-term fix is to reconfigure the machine (and
rebuild fesm-utils) onto a single netCDF-Fortran ≥ 4.6 for everything.

Note that the requirement is a property of *the archive you install*, not of
VILMA in general — an older VILMA build recipe found on this machine
(`/home/m/m300792/ModelCodes/VILMA3D/Makefile`) compiles against 4.5.3, but the
`vega_pism.a` actually installed at `VILMAROOT` here was built elsewhere against
a newer netCDF. If you install a different VILMA, re-check with the `dv.f90`
reproducer above before assuming either version.

If the link ever fails with `relocation truncated to fit: R_X86_64_32S`, add
`-mcmodel=large -shared-intel` (VILMA's own build uses them); the current install
links cleanly without them.

### 6b. VILMA must start from the reference state

VILMA requires the first record of its load history to vanish against the
reference ice file, and aborts otherwise:

```
 error: first referenced load file should vanish but shows range of
  -632.010925640728        725.757931451403
```

So `fe_vilma` pins history slice 1 to the reference ice (`h_ice_eq`) and never
rewrites it, exactly as CLIMBER-X does. The consequence is real and worth
stating: **VILMA begins with zero load anomaly and zero viscous memory at `t0`.**
If the run's start-slice ice is not the reference ice — the normal case for
`i_eq=1`, a present-day reference with an LGM start — VILMA absorbs the whole
start-vs-reference difference as a jump in the first interval, while the native
solver measures a genuine departure from a relaxed reference. `fe_vilma_init`
prints a WARNING when the two differ by more than 1 m.

For a like-for-like comparison, make the reference *be* the start state: `&ctl`
`i_eq=0` where the forcing file carries a bed for the start slice, or `i_eq=2`
with `h_ice_eq_file` set to the start-slice ice.

### 6c. Load ramping within an interval differs

With `vg%l_load_hist = .false.` the ice-history file keeps two slices and VILMA
re-reads them every step. Slice 1 is pinned to the reference (6b) and the epoch
coordinate is written once, so VILMA's internal time runs past the labelled
window after the first interval and it then holds the load at the current slice.
The native solver instead **ramps the load linearly** from the previous slice to
the new one across the interval. Same end points, different path within a
coupling step; expect first-order-in-Δt differences from this alone.

### 6d. Other operational notes

* **VILMA is a singleton.** All of its state is module-level `SAVE`d data, so
  only one backend instance per process, and it is not thread-safe at the API
  level.
* **Fixed unit numbers.** VILMA hardwires Fortran units 16, 23, 44–47, 55–58,
  66, 71–78, 80–93, 99. Do not open files on those from the host.
* **Relative default paths.** Every `io_*%n` defaults to a relative `inp/`,
  `out/` or `restart/` path. `fe_vilma` sets them all explicitly — including
  `io_rsrpt` (the rotation-potential state), which CLIMBER-X's wrapper forgets,
  leaving it writing to `./restart/rotpot.log` relative to the process cwd.
  `io_lis` is the one exception that still lands in the cwd as `vega.lis`.
* **120-character paths.** VILMA stores file names in `character(len=120)`.
  `fe_vilma_init` checks `vilma_out_dir` and fails loudly rather than letting a
  long path be truncated. Prefer a short relative `vilma_out_dir`.
* **`setup` is deferred.** It needs the time step and the start epoch, neither of
  which the `fe_coupling` API supplies at init, so `fe_vilma` calls it on the
  first *advancing* update and takes `vg%dt` from that interval. A forcing axis
  with a non-uniform cadence would therefore hand VILMA a step it was not
  configured for.
* **`openmp=0` will not link** with `vilma=1`: VILMA needs the OpenMP runtime
  (`__kmpc_fork_call`, `omp_get_wtime`). Build with `openmp=1` (the default).
* **VILMA's clock** is fed the model time in kyr (the forcing axis / 1000). Its
  response depends only on time differences, so the labelling affects only its
  own diagnostic output files.
* **`SLI_data.inp`** has no code references in this build; it is inert.

---

## 7. Verified run

Low-resolution smoke test, `lmax=32`, `jmax=170`, GLAC-1D forcing, `-1000 → 0` yr
at 100 yr (10 coupling steps), 1-D viscosity, on a levante login node with 8
threads:

```
   t=     -100.00 yr   max|rsl|=  2.98E+01
   t=        0.00 yr   max|rsl|=  3.17E+01
 [PROFILE] per coupling step (mean over 10 steps):
   read_ice (remap+IO) =    10.5 ms (  1.5 %)
   solid_earth_update  =   680.1 ms ( 97.4 %)
   fe_write_step (out) =     7.8 ms (  1.1 %)

 [PROFILE] solid_earth_update breakdown (per step, wall-clock):
   VILMA time_evolution  =   675.8 ms ( 99.4 % of update)
   Gauss<->VILMA remap   =     4.1 ms (  0.6 % of update)
   other (diagnostics)   =     0.1 ms (  0.0 % of update)
 fastearth: wrote out.nc
```

`out.nc` carries the usual `h_ice`, `rsl`, `z_bed`, `C_ocean`, `bsl` on the
128 × 66 model Gauss grid, 11 time slices; `rsl` grows monotonically from 0 to
−31.7 … +20.2 m and `bsl` stays near +2.2 m.
