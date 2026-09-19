# Reference state (present-day topography + ice) on the Gauss grid

Canonical present-day reference fields for the driver's `i_eq=1` (and `=3`) modes.
With `i_eq=1` the relaxed reference `z_bed_eq` (= the SLE `topo0`) and `h_ice_eq`
are this present-day state, so relative sea level is measured against today and
`rsl ≈ 0` at the present day.

There is a single canonical reference, **`rtopo_gauss_l128.nc`** (512×258 Gauss
grid). Point every run at it regardless of resolution:

```
z_bed_ref_file = "data/reference/rtopo_gauss_l128.nc"
h_ice_ref_file = "data/reference/rtopo_gauss_l128.nc"
```

`read_ref2d` (fe_drive) conservatively remaps it onto the run's own Gauss grid when
the resolutions differ, and reads it directly at lmax=128. The remap weights are
cached on disk (`fe_remap` → coords, under `maps/`), so the build cost is paid once
and reused on later runs and restarts — there is no longer any need for
per-resolution reference files. Because all resolutions derive from the same source,
their inputs stay directly comparable. (lmax≳128 is past the real structure in the
data, so l128 is the practical ceiling.)

## Provenance

Source: `RTopo-2.0.1_0.125deg_DRThydrocorr.nc` (CLIMBER-X input; `bedrock_topography`,
`ice_thickness`). Offline pipeline:

1. `rtopo_0.5deg.nc` — the 0.125° RTopo coarsened to 0.5° by an exact 4×4 block mean
   (conservative on a regular grid; the 0.5° `topo_05x05.nc` in CLIMBER-X carries bed
   only, so the full RTopo is coarsened here to keep ice too). Kept as the master
   source.
2. `rtopo_gauss_l128.nc` — generated from the 0.5° source by `fastearth_mkref`
   (`bin/fastearth_mkref.x mkref_l128.nml fastearth.nml`), conservatively remapping
   bed (as-is) and ice (mass-conserving) with the same `fe_remap` engine the online
   path uses.

To regenerate (or produce a different ceiling resolution), edit `mkref_l128.nml` and
run `make fastearth_mkref` then `bin/fastearth_mkref.x mkref_l128.nml fastearth.nml`.
