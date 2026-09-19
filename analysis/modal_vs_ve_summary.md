# Modal vs VE — accuracy & cost summary

Each modal run is compared against the VE run in the same set. `rsl RMSE` is area-weighted (cos lat) over all space-time; `rsl PD` at present day; `cost` is the per-coupling-step solver time; `speedup` = VE/modal.

### radial

VE reference: solver cost = 61.2 ms/step, n_solve = 1.6 sub-steps/step

VE scales: |rsl|max(PD) = 146.0 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 0.812 m | 0.342 m | 7.790 m | 0.009 m | 0.000 m | 53.2 ms | 1.2x | 10.3 |
| n1/isostatic | 17.453 m | 6.049 m | 195.197 m | 0.102 m | 0.002 m | 34.1 ms | 1.8x | 5.2 |
| n1/rate | 19.354 m | 6.196 m | 213.002 m | 0.153 m | 0.002 m | 55.4 ms | 1.1x | 8.7 |
| n1/residue | 23.962 m | 7.217 m | 271.688 m | 0.127 m | 0.004 m | 22.9 ms | 2.7x | 4.5 |
| n2/isostatic | 6.081 m | 3.810 m | 61.558 m | 0.041 m | 0.002 m | 37.2 ms | 1.6x | 7.5 |
| n2/rate | 5.555 m | 4.927 m | 51.273 m | 0.039 m | 0.001 m | 36.3 ms | 1.7x | 7.3 |
| n2/residue | 7.657 m | 4.957 m | 62.168 m | 0.074 m | 0.001 m | 40.4 ms | 1.5x | 8.1 |
| n4/isostatic | 0.919 m | 0.360 m | 7.872 m | 0.011 m | 0.000 m | 41.8 ms | 1.5x | 7.9 |
| n4/rate | 1.043 m | 0.723 m | 8.287 m | 0.012 m | 0.005 m | 49.8 ms | 1.2x | 9.9 |
| n4/residue | 0.933 m | 0.387 m | 7.865 m | 0.011 m | 0.000 m | 50.4 ms | 1.2x | 9.6 |
| n8/isostatic | 0.813 m | 0.353 m | 7.764 m | 0.009 m | 0.000 m | 46.0 ms | 1.3x | 9.0 |
| n8/rate | 0.812 m | 0.341 m | 7.808 m | 0.009 m | 0.000 m | 54.2 ms | 1.1x | 10.8 |
| n8/residue | 0.813 m | 0.342 m | 7.793 m | 0.009 m | 0.000 m | 53.3 ms | 1.1x | 10.6 |

### deglac3d

VE reference: solver cost = 1756.3 ms/step, n_solve = 8.6 sub-steps/step

VE scales: |rsl|max(PD) = 342.6 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 21.053 m | 23.513 m | 402.553 m | 0.188 m | 0.000 m | 1359.3 ms | 1.3x | 19.2 |
| n1/isostatic | 32.471 m | 28.246 m | 509.605 m | 0.242 m | 0.002 m | 218.2 ms | 8.0x | 20.3 |
| n1/rate | 33.398 m | 28.243 m | 510.927 m | 0.273 m | 0.002 m | 338.3 ms | 5.2x | 31.9 |
| n1/residue | 36.448 m | 28.344 m | 522.573 m | 0.259 m | 0.006 m | 187.1 ms | 9.4x | 17.2 |
| n2/isostatic | 24.758 m | 26.586 m | 454.866 m | 0.208 m | 0.001 m | 236.3 ms | 7.4x | 14.4 |
| n2/rate | 25.663 m | 28.146 m | 471.229 m | 0.207 m | 0.007 m | 326.0 ms | 5.4x | 19.8 |
| n2/residue | 27.056 m | 28.258 m | 474.792 m | 0.230 m | 0.007 m | 318.8 ms | 5.5x | 19.4 |
| n4/isostatic | 21.142 m | 23.547 m | 403.918 m | 0.190 m | 0.000 m | 394.5 ms | 4.5x | 14.3 |
| n4/rate | 21.798 m | 24.258 m | 418.568 m | 0.190 m | 0.001 m | 537.4 ms | 3.3x | 19.9 |
| n4/residue | 22.424 m | 24.204 m | 419.983 m | 0.193 m | 0.002 m | 357.9 ms | 4.9x | 13.3 |
| n8/isostatic | 21.061 m | 23.521 m | 402.735 m | 0.188 m | 0.000 m | 936.3 ms | 1.9x | 19.2 |
| n8/rate | 21.653 m | 24.040 m | 416.703 m | 0.187 m | 0.000 m | 910.8 ms | 1.9x | 19.1 |
| n8/residue | 22.308 m | 24.142 m | 418.699 m | 0.191 m | 0.003 m | 1252.0 ms | 1.4x | 25.5 |

