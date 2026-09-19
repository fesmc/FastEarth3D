#!/usr/bin/env julia
# Analysis plots for the §14 forced deglaciation runs (LGM -> PD, Tarasov + RTopo
# PD reference, i_eq=1). For each run: a 3-panel map (-20 ka, -10 ka, PD) showing
#   - bedrock topography (semi-transparent topographic colors),
#   - ice-sheet surface shaded grey->white with 500 m surface contours on top,
#   - relative-sea-level contours (labelled) over the ocean,
#   - the barystatic sea level (BSL) printed in the top-right corner.
# Plus one figure of the BSL time series for every run on the same axes.
#
# Usage:  julia analysis/plot_deglac.jl [run1=out1.nc run2=out2.nc ...]
# With no args it uses the default RUNS table below.

using NCDatasets
using CairoMakie
using Printf

# ---------------------------------------------------------------------------
# runs: label => output netCDF. Override on the command line as label=path.
const DEFAULT_RUNS = [
    "Bagge2021"      => "runs/s14_bagge/out.nc",
    "Pan2022"        => "runs/s14_pan/out.nc",
    "Bagge +1σ"      => "runs/s14_bagge_p1/out.nc",
    "Bagge -1σ"      => "runs/s14_bagge_m1/out.nc",
]

const TIMES_KA = [-20.0, -10.0, 0.0]      # panel times [ka]
const OUTDIR   = "analysis/figs"

# topographic colormap: blues below sea level, greens/browns/white above.
const TOPO_CMAP = cgrad([:navy, :dodgerblue, :paleturquoise,
                         :darkgreen, :khaki, :saddlebrown, :white],
                        [0.0, 0.30, 0.48, 0.52, 0.68, 0.88, 1.0])
const TOPO_RANGE = (-6000.0, 6000.0)
const ICE_CMAP   = cgrad([:gray35, :gray70, :white])
const ICE_RANGE  = (0.0, 4000.0)

# ---------------------------------------------------------------------------
parse_runs(args) = isempty(args) ? DEFAULT_RUNS :
    [ (p = split(a, "="; limit = 2); String(p[1]) => String(p[2])) for a in args ]

"Roll lon to [-180,180), sort lon and lat ascending; return (lon, lat, perm_lon, perm_lat)."
function grid_axes(ds)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    lon = map(x -> x > 180 ? x - 360 : x, lon)
    pl = sortperm(lon); pt = sortperm(lat)
    return lon[pl], lat[pt], pl, pt
end

"2-D field name at time index it, reordered to ascending lon/lat."
field(ds, name, it, pl, pt) = Float64.(ds[name][:, :, it])[pl, pt]

nearest_index(times_yr, t_ka) = argmin(abs.(times_yr ./ 1000 .- t_ka))

function panel!(ax, ds, it, pl, pt, lon, lat)
    zb  = field(ds, "z_bed",   it, pl, pt)
    hi  = field(ds, "h_ice",   it, pl, pt)
    C   = field(ds, "C_ocean", it, pl, pt)
    surf = zb .+ hi                                   # ice-surface elevation
    has_ice = hi .> 1.0

    # base bedrock topography (semi-transparent)
    heatmap!(ax, lon, lat, zb; colormap = TOPO_CMAP, colorrange = TOPO_RANGE,
             alpha = 0.65)

    # ice surface shaded grey->white where grounded ice exists
    ice_s = map((s, m) -> m ? s : NaN, surf, has_ice)
    heatmap!(ax, lon, lat, ice_s; colormap = ICE_CMAP, colorrange = ICE_RANGE)

    # 500 m ice-surface contours on top of the ice
    contour!(ax, lon, lat, ice_s; levels = 500:500:4500,
             color = (:black, 0.45), linewidth = 0.5)

    # RSL contours over the ocean only (labelled)
    is_sea = (C .> 0.5) .| (zb .< 0.0)
    rsl = field(ds, "rsl", it, pl, pt)
    rsl_sea = map((r, m) -> m ? r : NaN, rsl, is_sea)
    contour!(ax, lon, lat, rsl_sea; levels = -160:40:160, labels = true,
             labelsize = 8, color = :firebrick, linewidth = 0.8)

    # BSL in the top-right corner
    bsl = Float64(ds["bsl"][it])
    text!(ax, 0.98, 0.97; text = @sprintf("BSL = %.0f m", bsl),
          space = :relative, align = (:right, :top), fontsize = 11,
          font = :bold, color = :black,
          strokecolor = :white, strokewidth = 2)
    return nothing
end

function plot_run(label, path)
    ds = NCDataset(path)
    times = Float64.(ds["time"][:])              # years
    lon, lat, pl, pt = grid_axes(ds)

    fig = Figure(size = (1500, 520))
    Label(fig[0, 1:3], label; fontsize = 18, font = :bold)
    for (j, tka) in enumerate(TIMES_KA)
        it = nearest_index(times, tka)
        ax = Axis(fig[1, j]; aspect = DataAspect(),
                  title = @sprintf("%.0f ka  (t = %.0f yr)", tka, times[it]),
                  xlabel = "lon", ylabel = j == 1 ? "lat" : "")
        limits!(ax, -180, 180, -90, 90)
        panel!(ax, ds, it, pl, pt, lon, lat)
    end
    Colorbar(fig[1, 4], colormap = TOPO_CMAP, limits = TOPO_RANGE,
             label = "bedrock elevation [m]")
    close(ds)
    fn = joinpath(OUTDIR, "maps_" * replace(label, r"[^A-Za-z0-9]" => "_") * ".png")
    save(fn, fig); println("wrote ", fn)
end

function plot_bsl(runs)
    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time [ka]", ylabel = "barystatic sea level [m]",
              title = "Barystatic sea level vs present-day reference")
    for (label, path) in runs
        isfile(path) || continue
        ds = NCDataset(path)
        t = Float64.(ds["time"][:]) ./ 1000
        b = Float64.(ds["bsl"][:])
        # drop the t0 seed sample (reference initialisation, pre-SLE-solve)
        lines!(ax, t[2:end], b[2:end]; label = label, linewidth = 2)
        close(ds)
    end
    axislegend(ax; position = :rb)
    fn = joinpath(OUTDIR, "bsl_timeseries.png")
    save(fn, fig); println("wrote ", fn)
end

function main()
    runs = parse_runs(ARGS)
    mkpath(OUTDIR)
    for (label, path) in runs
        isfile(path) || (@warn "missing output, skipping" label path; continue)
        plot_run(label, path)
    end
    plot_bsl(runs)
end

main()
