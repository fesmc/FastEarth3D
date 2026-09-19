# FastEarth3D

A state-of-the-art but simple and fast **3D solid-Earth model** — a
visco-elastic deformation model coupled with the sea-level equation — intended
as an **open-source replacement for VILMA** within the CLIMBER-X climate model.

The method is a clean-room reimplementation of the **spectral–finite-element,
time-domain** approach of Martinec (2000): spherical harmonics horizontally,
finite elements radially, an incompressible Maxwell rheology integrated
explicitly in time, a self-consistent sea-level equation with migrating
coastlines, and rotational feedback. It is built **3D-ready from the start**
(laterally varying viscosity) and validated against the published GIA benchmarks.

Full documentation (physics, discretization, benchmarks, install & run) is the
Quarto site under [`docs/`](docs/); see [doc/design.md](doc/design.md) for the
design rationale and method comparison.

## Status

The full model is implemented and validated: the spectral–finite-element solver
core, viscoelastic time stepping, the self-consistent migrating-coastline
sea-level equation, rotational feedback (polar motion), and laterally varying
(3D) viscosity — plus restart, spin-up, online lon-lat→Gauss remapping, and a
host-coupling API. A second response operator, a reduced **modal solver**, offers
a tunable, faster approximation that converges back to the full solver.

Validated against the Spada et al. (2011) and Martinec et al. (2018) community
benchmarks (radial Love numbers, disc-load response, sea-level equation) and
Spada test 3/2 (rotation). Cross-code validation of the 3D path and the
modal-vs-VE accuracy/cost study are ongoing.

## Install

Dependencies are [fesm-utils](https://github.com/fesm-org/fesm-utils) (branch
`coords-dev`, providing FFTW, SHTns, the `fesmutils` helper and the `coords`
module) plus a system netCDF. **configme** clones/links the dependencies and
generates the machine/compiler Makefile:

```bash
configme install FastEarth3D                 # resolve machine/compiler, clone deps
configme install FastEarth3D -m macbook -c gfortran
configme install FastEarth3D --link fesm-utils=/abs/path/to/fesm-utils   # reuse a checkout
```

## Build

```bash
make fastearth        # -> bin/fastearth.x        (standalone forced-run driver)
make fastearth_mkref  # -> bin/fastearth_mkref.x  (build a Gauss-grid reference)
make fastearth_remap  # -> bin/fastearth_remap.x  (offline lon-lat -> Gauss remap)
make check            # build + run the test suite
```

`make` switches: `debug=0|1|2`, `openmp=0|1` (default 1; build `openmp=1` for the
threaded degree loop at production resolutions).

## Configure & run

All runtime parameters live in a single namelist group `&fe3d`, loaded into the
`fe_param_class` record by `fe_par_load`. [`fastearth.nml`](fastearth.nml) is the
complete, documented defaults set; a run can pass a sparse file overlaid on it
(yelmo `defaults_file` convention), overriding only what it needs. Time fields
(`dt_*`, `time_*`) are given in **years** and converted to SI seconds on load.

- **Earth structure** — `earth`: a named built-in (e.g. `"M3-L70-V01"`) or
  `"custom"` to assemble from the surface-first layer arrays.
- **Response solver** — `earth_response`: `"ve"` (full viscoelastic, default),
  `"modal"` (reduced, with `n_modes` / `mode_rank`), `"elastic"`, `"null"`.
- **Time scheme** — `scheme = "fe"` (1st-order explicit) or `"trap"` (2nd-order
  adaptive), advanced by the `fe_timestep` controller.
- **3D viscosity / spin-up / restart** — `l_visc_3d`, `dt_equil`, `spinup_1d`,
  `restart_in_file`.
- **Rotation** — `rotation` (TPW feedback): on by default; `.false.` for the
  non-rotating benchmarks.

Run the standalone driver directly (sparse overlay + complete defaults):

```bash
./bin/fastearth.x examples/deglac_lgm.nml fastearth.nml
```

It reads a reference state and an ice-thickness forcing (`file_forcing`,
`h_ice(lon,lat,time)`; remapped from lon-lat on the fly by default), marches the
model across the forcing, and writes the diagnostic surface fields (`rsl`,
`z_bed`, …) to `file_out`.

Or stage/submit runs and ensembles with **runme** (`-r` run, `-s` submit;
comma-lists in `-p` define ensemble dimensions):

```bash
runme -o runs/deglac -e main --omp 8 -r -p fe3d.lmax=128 fe3d.earth_response=ve
```

The [`scripts/run_modal_vs_ve.sh`](scripts/run_modal_vs_ve.sh) launcher stages
the full modal-vs-VE accuracy/cost sweep through runme.

Embedding the model in a host (the CLIMBER-X coupling path) uses the same API
behind a single `use fastearth3d`:

```fortran
use fastearth3d
type(fe_param_class) :: par
type(solid_earth)    :: se
call fe_par_load(par, "fastearth.nml")
call solid_earth_init(se, par, sht, z_bed_eq, h_ice_ref)
call solid_earth_update(se, h_ice, dt)   ! advance time -> time+dt; reads se%rsl, se%z_bed
```
