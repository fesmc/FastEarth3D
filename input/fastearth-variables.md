# fastearth restart / output variables

Variable-io table (yelmo convention) for FastEarth3D netCDF I/O, used by
`fe_io` for both restart files and diagnostic `write_step` output. Each row gives
the netCDF variable name, its dimensions, units, and a long_name. The time axis
(unlimited) lets several snapshots live in one file.

The prognostic state restored on restart depends on the response kind: the
Maxwell memory-stress fields `tau_*` (+ `sigma_n_*`) for `RESP_VE`, or the modal
amplitudes `phi_*` for `RESP_MODAL`; plus the adaptive controller's next-step Δt
seed `dt_try` for both. The reference fields `z_bed_eq`/`h_ice_eq` are static
(written once) and checked on read; the rest are diagnostic.

| id | variable     | dimensions    | units  | long_name                                        |
|----|--------------|---------------|--------|--------------------------------------------------|
|  1 | tau_a_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (real part)   |
|  2 | tau_a_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (imag part)   |
|  3 | tau_b_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (real part)   |
|  4 | tau_b_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (imag part)   |
|  5 | tau_c_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (real part)   |
|  6 | tau_c_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (imag part)   |
|  7 | z_bed_eq     | lon, lat      | m      | Reference (equilibrium) bedrock elevation        |
|  8 | h_ice_eq     | lon, lat      | m      | Reference (equilibrium) grounded-ice thickness   |
|  9 | h_ice        | lon, lat      | m      | Grounded-ice thickness                           |
| 10 | rsl          | lon, lat      | m      | Relative sea-level change (full field)           |
| 11 | z_bed        | lon, lat      | m      | Bedrock elevation (z_bed_eq - rsl)               |
| 12 | C_ocean      | lon, lat      | 1      | Ocean function (1 ocean / 0 land)                |
| 13 | dt_try       | time          | years  | Adaptive time-stepping next-step Δt suggestion   |
| 14 | sigma_n_re   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (real part)   |
| 15 | sigma_n_im   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (imag part)   |
| 16 | sigma_primed | time          | 1      | Flag: σ_n is tracked (1) or not yet primed (0)   |
| 17 | bsl          | time          | m      | Barystatic sea level vs reference (eustatic eq.) |
| 18 | phi_re       | nphi_modal, time | Pa  | Modal amplitude φ per (l,m) and mode (real part)  |
| 19 | phi_im       | nphi_modal, time | Pa  | Modal amplitude φ per (l,m) and mode (imag part)  |
| 20 | rot_m_re     | time          | 1      | Polar motion m₁ (rotation, real part)            |
| 21 | rot_m_im     | time          | 1      | Polar motion m₂ (rotation, imag part)            |
| 22 | rot_load_mem | nlam, ne_rot, nrc | Pa | Rotation loading-channel memory stress (packed)  |
| 23 | rot_tidal_mem| nlam, ne_rot, nrc | Pa | Rotation tidal-channel memory stress (packed)    |
