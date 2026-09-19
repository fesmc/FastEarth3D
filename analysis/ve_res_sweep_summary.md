# VE resolution & sub-step sweep — accuracy & cost

Reference: lmax=128, cfl=1. Every run is regridded onto the reference grid; `rsl RMSE` is area-weighted (cos lat) over all space-time, `rsl PD` at present day.

Reference cost: solid_earth_update = 13321.0 ms/step (drift 357.0 + memory 12108.6), n_solve = 11.9

Reference scales: |rsl|max(PD) = 289.2 m, BSL min = -120.1 m, BSL(PD) = 1.44 m

| run | lmax | cfl | vtol | ne3d | rsl RMSE | rsl PD | rsl max | bsl PD | cost | mem | nsolve | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lmax.32 | 32 | 1 | 1e-03 | 51 | 11.242 m | 5.129 m | 168.670 m | 0.757 m | 327.0 ms | 222.0 | 11.9 | 40.74x |
| lmax.64 | 64 | 1 | 1e-03 | 51 | 2.918 m | 2.365 m | 59.523 m | 0.664 m | 1763.7 ms | 1281.8 | 11.9 | 7.55x |
| lmax.96 | 96 | 1 | 1e-03 | 51 | 1.571 m | 1.544 m | 38.069 m | 0.556 m | 4165.4 ms | 3430.0 | 11.9 | 3.20x |
| lmax.128 | 128 | 1 | 1e-03 | 66 | 0.000 m | 0.000 m | 0.000 m | 0.000 m | 13321.0 ms | 12108.6 | 11.9 | 1.00x |
| lmax.128.cfl0.5 | 128 | 0.5 | 1e-03 | 66 | 0.024 m | 0.015 m | 0.528 m | 0.000 m | 25858.4 ms | 23523.2 | 22.8 | 0.52x |
| lmax.128.cfl1.5 | 128 | 1.5 | 1e-03 | 66 | 0.045 m | 0.021 m | 0.837 m | 0.000 m | 8294.5 ms | 7517.7 | 7.6 | 1.61x |
| lmax.128.cfl2.0 | 128 | 2 | 1e-03 | 66 | 0.058 m | 0.028 m | 1.294 m | 0.000 m | 6669.3 ms | 6042.5 | 6.0 | 2.00x |
| lmax.128.cfl2.5 | 128 | 2.5 | 1e-03 | 66 | NaN m | NaN m | NaN m | 1.436 m | 5018.2 ms | 4546.3 | 13.1 | 2.65x |
| lmax.128.vtol0.1 | 128 | 1 | 1e-01 | 65 | 0.001 m | 0.001 m | 0.040 m | 0.000 m | 13100.4 ms | 11883.5 | 11.9 | 1.02x |
| lmax.128.vtol0.3 | 128 | 1 | 3e-01 | 62 | 0.029 m | 0.028 m | 0.900 m | 0.000 m | 12430.4 ms | 11212.3 | 11.9 | 1.07x |
| lmax.128.vtol1.0 | 128 | 1 | 1e+00 | 52 | 0.078 m | 0.074 m | 2.391 m | 0.000 m | 10717.7 ms | 9501.4 | 11.9 | 1.24x |

