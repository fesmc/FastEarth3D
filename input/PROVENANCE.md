# `input/` data provenance

## `geo_ice_glac1d_deglac.nc`

Ice-thickness forcing for the **last-deglaciation experiment**: GLAC-1D
(Tarasov) ice thickness, LGM → present-day.

- **Source:** `geo_ice_tarasov_deglac.nc` as distributed with **CLIMBER-X**
  (`input/`), labelled "Tarasov PMIP4" per-variable and "Tarasov GLAC1D" in the
  companion `geo_ice_tarasov_lgc_0ka.nc`.
- **Underlying dataset:** GLAC-1D (Tarasov, Briggs & Peltier), the PMIP4
  deglaciation ice reconstruction.
- **Grid:** 360 lon × 360 lat (1° × 0.5°), 261 time slices, **−26 000 → 0 yr**
  at 100 yr. That is the full span of the source file — no time subsetting was
  applied.
- **Contents:** `ice_thickness` only.

### Why only ice thickness, and how it got to 7 MB

The source file is **1.08 GB**; this one is **7.1 MB**, a factor of ~152. Three
independent reductions, none of which touch what the model reads:

1. **One variable.** The standalone driver's forcing is `ice_thickness`
   (`fe3d.name_ice`); the present-day reference topography comes from RTopo-2
   (`data/reference/rtopo_gauss_l128.nc`), not from this file.
   `surface_elevation` is bedrock+ice and `mask` is derivable, so neither is
   independent data. `bedrock_topography` is the expensive one — 35–45 MB even
   packed, against 8.8 MB for the ice — and the model does not read it here. It
   remains in the full CLIMBER-X file if a comparison ever needs it.
2. **int16 at 1 m.** Ice thickness spans 0–4108 m, so 1 m resolution fits int16
   with wide margin and is far finer than the reconstruction's own accuracy.
   Maximum round-off **0.500 m**; total ice volume is preserved to **0.0001 %**.
3. **Deflate (level 6) + shuffle.** The field is 84 % zeros (ice-free), which
   compresses hard. Chunking across time was measured and does **not** help —
   deflate does not exploit the slab-to-slab redundancy — so chunks stay at
   4 time steps.

The source time axis carries float noise (`-899.999999999991`); it is rounded to
whole years here.

### Regenerating

Built by `experiments/make_glac1d_forcing.py` in the
`paper-fastearth3d-experiments` repository, which also records the exact command
and the output checksum in its `LOG.md`. The script takes `--t0/--t1` to narrow
the window (e.g. `--t0 -21000` for the 21 ka LGM-snapshot convention) and
`--with-bedrock` to include the bedrock field.

```
sha256  65a10e3acfcaa6219d4817494c91d48a174a74df6548184f927233030597d315
```

## `bagge2021.nc`, `pan2022.nc`

Laterally varying (3-D) mantle viscosity fields, `log10(eta)`. These predate
this file and their upstream provenance is not recorded here — it should be
filled in by whoever staged them.
