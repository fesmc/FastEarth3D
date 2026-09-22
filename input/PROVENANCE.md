# `input/` data provenance

## `geo_ice_glac1d_deglac.nc`

Forcing for the **last-deglaciation experiment**: GLAC-1D
(Tarasov) ice thickness, LGM → present-day, with bedrock for comparison.

- **Source:** `geo_ice_tarasov_deglac.nc` as distributed with **CLIMBER-X**
  (`input/`), labelled "Tarasov PMIP4" per-variable and "Tarasov GLAC1D" in the
  companion `geo_ice_tarasov_lgc_0ka.nc`.
- **Underlying dataset:** GLAC-1D (Tarasov, Briggs & Peltier), the PMIP4
  deglaciation ice reconstruction.
- **Grid:** 360 lon × 360 lat (1° × 0.5°), 261 time slices, **−26 000 → 0 yr**
  at 100 yr. That is the full span of the source file — no time subsetting was
  applied.
- **Contents:** `ice_thickness`, plus `bedrock_topography_pd` (present-day
  bedrock, float64) and `bedrock_anomaly` (bedrock change relative to it).

### Bedrock is stored as a change relative to present day

```
bedrock(t) = bedrock_topography_pd + bedrock_anomaly(t)
```

The model reads **neither** — its reference topography is RTopo-2
(`data/reference/rtopo_gauss_l128.nc`). Bedrock is here so simulated bedrock
change can be compared against the reconstruction, and the anomaly is the
*directly comparable* quantity: a GIA model predicts change relative to a
present-day reference, which is what `z_bed` and `rsl` in the model output are.

It is also much smaller. Absolute bedrock packs to **45.4 MB**; the anomaly to
**11.6 MB**, a factor of 3.9 at identical precision and cadence. The static
bathymetry detail that dominates the entropy is byte-identical in every slice
and cancels, leaving only the smooth GIA signal. `bedrock_topography_pd` is kept
at full float64 (it costs ~0.3 MB), so the reconstruction above is exact to the
anomaly's own 0.5 m rounding, and the present-day slice is bit-exact.

### How it got from 1.08 GB to 19 MB

A factor of ~58. None of these reductions touch what the model reads:

1. **Only what is needed.** `surface_elevation` is bedrock+ice and `mask` is
   derivable, so neither is independent data; both are dropped.
2. **int16 at 1 m.** Ice thickness spans 0–4108 m, so 1 m resolution fits int16
   with wide margin and is far finer than the reconstruction's own accuracy.
   Maximum round-off **0.500 m** for both ice and the bedrock anomaly; total ice
  volume is preserved to **0.0001 %**.
3. **Deflate (level 6) + shuffle.** Ice is 84 % zeros (ice-free), which
   compresses hard. Chunking across time was measured and does **not** help —
   deflate does not exploit the slab-to-slab redundancy — so chunks stay at
   4 time steps.
4. **Bedrock as an anomaly**, as above.

Measured storage rates, useful for sizing any other window:

| Variable | Size | Per slice |
|---|---:|---:|
| `ice_thickness` | 7.1 MB | 27.4 kB |
| `bedrock_anomaly` | 11.6 MB | 44.5 kB |
| bedrock stored absolutely (rejected) | 45.4 MB | 174.0 kB |

A **full last glacial cycle** (122 kyr) at the same 100 yr cadence would be
**1221 slices**: ~33 MB for ice alone, ~88 MB with bedrock. At 500 yr it drops
to ~7 MB / ~18 MB. Those are extrapolations from the rates above — the 122 kyr
GLAC-1D series is **not** staged on Levante (`geo_ice_tarasov_lgc_0ka.nc` holds
a single present-day slice), and the glacial period has more ice cover and
larger bedrock anomalies than the deglaciation, so treat them as lower bounds.

The source time axis carries float noise (`-899.999999999991`); it is rounded to
whole years here.

### Regenerating

Built by `experiments/make_glac1d_forcing.py` in the
`paper-fastearth3d-experiments` repository, which also records the exact command
and the output checksum in its `LOG.md`. The script takes `--t0/--t1` to narrow
the window (e.g. `--t0 -21000` for the 21 ka LGM-snapshot convention) and
`--with-bedrock` to include the bedrock fields (used for the shipped file).

```
sha256  ed4c141605d914140648d116f61fa6e8b7351105037af193a7556fc46e7de72f
```

## `bagge2021.nc`, `pan2022.nc`

Laterally varying (3-D) mantle viscosity fields, `log10(eta)`. These predate
this file and their upstream provenance is not recorded here — it should be
filled in by whoever staged them.

`bagge2021.nc` stores latitude NORTH-first (descending). That is legitimate —
`fe_read_visc_3d` normalises either orientation to ascending on read — but note
it when comparing against a backend that reads the file itself.

`pan2022.nc` has been edited since staging: longitude index 360 duplicated
index 359 (both 0.0 deg, with a bit-identical `eta` column), which left the
axis non-monotonic and unusable for interpolation. The duplicate column was
dropped, giving 720 longitudes on a uniform 0.5 deg grid spanning exactly
360 deg; `lat`, `r` and every retained `eta` value are untouched. **The same
defect is present upstream** in
`isostasy_data/earth_structure/viscosity/pan2022.nc`, so re-staging from there
will reintroduce it.
