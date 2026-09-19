#!/usr/bin/env julia
# Spin-up comparison plots (branch restart-spinup). For each run: a 2x2 map grid
# at [init after spin-up, -20 ka, -10 ka, 0 ka] showing
#   - bedrock topography (semi-transparent topographic colors) on land/shelf only,
#   - ice-sheet surface shaded grey->white with 500 m surface contours (black),
#   - relative sea level as a shaded overlay + grey contours over the ocean,
#   - the barystatic sea level (BSL) printed plain black in the top-right corner.
# Three colorbars (bedrock, ice surface, relative sea level) sit on the right.
# Plus one figure of the BSL time series for every run on the same axes.
#
# Usage:  julia analysis/plot_spinup_comparison.jl [label=out.nc ...]
# With no args it uses the RUNS table below (the l64 1-D vs 3-D spin-up experiment
# plus the cross-resolution l128-from-l64 run).

using NCDatasets
using CairoMakie
using Printf

# label => output netCDF. Override on the command line as label=path.
const DEFAULT_RUNS = [
    "3-D spin-up (l64)"        => "runs/cmp_spin3d_l64/out.nc",
    "1-D spin-up (l64)"        => "runs/cmp_spin1d_l64/out.nc",
    "l128 from l64 spin-up"    => "runs/xres_l128_from_l64/out.nc",
    "1-D transient (l128)"     => "runs/trans1d_l128/out.nc",
]

# BSL time series: only these two runs, from BSL_TMIN_KA onward (skips the
# first-step artifact at -26 ka).
const BSL_RUNS = [
    "l128 from l64 spin-up"    => "runs/xres_l128_from_l64/out.nc",
    "1-D transient (l128)"     => "runs/trans1d_l128/out.nc",
]
const BSL_TMIN_KA = -25.9

# Four panels at the nearest slices to these times [ka]. Panel 1 uses -25.9 ka
# (not the -26 ka init slice, which has an artifact in one run).
const TIMES_KA = [-25.9, -20.0, -10.0, 0.0]
const OUTDIR   = "analysis/figs"

# Out-of-range triangles always take the palette's own end colours.
clip_low(cm)  = cm[0.0]
clip_high(cm) = cm[1.0]

# Bedrock colormap limited to [-1000, 4000] m with the ocean->land colour
# transition placed at 0 m (fraction 0.20 of the 5000 m span).
const TOPO_CMAP  = cgrad([:navy, :dodgerblue, :paleturquoise,
                          :darkgreen, :khaki, :saddlebrown, :white],
                         [0.0, 0.10, 0.18, 0.22, 0.50, 0.80, 1.0])
const TOPO_RANGE = (-1000.0, 4000.0)
const TOPO_LOW   = clip_low(TOPO_CMAP)
const TOPO_HIGH  = clip_high(TOPO_CMAP)
const TOPO_ALPHA = 0.5

const ICE_CMAP   = cgrad([:gray35, :gray70, :white])
const ICE_RANGE  = (0.0, 4000.0)
const ICE_HIGH   = clip_high(ICE_CMAP)

# Relative sea level: diverging overlay over the ocean with log-like colour steps.
# The labelled levels below are mapped to EVENLY spaced positions on the colour
# axis (piecewise-linear between breakpoints); data is plotted in that index
# space so each segment gets equal colour range. White is aligned with 0 m.
const RSL_TICKS  = Float64[-500, -100, -50, -10, -5, -1, 0, 1, 5, 10, 50, 100, 500]
const RSL_N      = length(RSL_TICKS)
const RSL_CRANGE = (0.0, Float64(RSL_N - 1))   # index space: breakpoints at 0..N-1

"Map an rsl value [m] to its evenly-spaced index position; linear within segments,
extrapolated past the ends so out-of-range data triggers the clip triangles."
function rsl_index(x)
    isnan(x) && return NaN
    B = RSL_TICKS; n = RSL_N
    x <= B[1] && return (x - B[1]) / (B[2] - B[1])
    x >= B[n] && return (n - 1) + (x - B[n]) / (B[n] - B[n-1])
    i = searchsortedlast(B, x)
    return (i - 1) + (x - B[i]) / (B[i+1] - B[i])
end

const _VIK       = cgrad(:vik)
const _RSL_WHITE = rsl_index(0.0) / (RSL_N - 1)
const RSL_CMAP   = cgrad([_VIK[0.0], RGBf(1, 1, 1), _VIK[1.0]], [0.0, _RSL_WHITE, 1.0])
const RSL_LOW    = clip_low(RSL_CMAP)
const RSL_HIGH   = clip_high(RSL_CMAP)
const RSL_CBTICKS = (collect(0.0:RSL_N-1), [@sprintf("%d", round(Int, b)) for b in RSL_TICKS])
const RSL_LEVELS  = Float64[-160, -140, -120, -100, -40, -20, -10, -5, 0, 5, 10]
const RSL_ALPHA   = 0.8

parse_runs(args) = isempty(args) ? DEFAULT_RUNS :
    [ (p = split(a, "="; limit = 2); String(p[1]) => String(p[2])) for a in args ]

"Roll lon to [-180,180), sort lon and lat ascending; return (lon, lat, perm_lon, perm_lat).
lon/lat are returned as evenly spaced ranges (required for heatmap interpolate=true).
The l128 Gauss latitudes are near-uniform in the interior (~0.696 deg), so the regular
axis only warps the mapping by a fraction of a degree near the poles."
function grid_axes(ds)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    lon = map(x -> x > 180 ? x - 360 : x, lon)
    pl = sortperm(lon); pt = sortperm(lat)
    slon = lon[pl]; slat = lat[pt]
    lonr = range(first(slon), last(slon), length = length(slon))
    latr = range(first(slat), last(slat), length = length(slat))
    return lonr, latr, pl, pt
end

field(ds, name, it, pl, pt) = Float64.(ds[name][:, :, it])[pl, pt]
nearest_index(times_yr, t_ka) = argmin(abs.(times_yr ./ 1000 .- t_ka))

function panel!(ax, ds, it, pl, pt, lon, lat)
    zb  = field(ds, "z_bed",   it, pl, pt)
    hi  = field(ds, "h_ice",   it, pl, pt)
    C   = field(ds, "C_ocean", it, pl, pt)
    surf = zb .+ hi
    has_ice  = hi .> 1.0
    is_ocean = C .> 0.5

    # Bedrock topography on land/exposed shelf only (no colours in the ocean).
    zb_land = map((z, o) -> o ? NaN : z, zb, is_ocean)
    heatmap!(ax, lon, lat, zb_land; colormap = TOPO_CMAP, colorrange = TOPO_RANGE,
             lowclip = TOPO_LOW, highclip = TOPO_HIGH, alpha = TOPO_ALPHA, interpolate = true)

    # Relative sea level: shaded overlay (log-like even-step colours) + grey
    # contours, ocean only. Heatmap is plotted in the even-spaced index space;
    # contours stay on the raw field so their labels read in metres.
    rsl = field(ds, "rsl", it, pl, pt)
    rsl_sea = map((r, o) -> o ? r : NaN, rsl, is_ocean)
    heatmap!(ax, lon, lat, rsl_index.(rsl_sea); colormap = RSL_CMAP, colorrange = RSL_CRANGE,
             lowclip = RSL_LOW, highclip = RSL_HIGH, alpha = RSL_ALPHA, interpolate = true)
    contour!(ax, lon, lat, rsl_sea; levels = RSL_LEVELS, labels = true,
             labelsize = 8, color = (:gray25, 0.8), linewidth = 0.6)

    # Ice-sheet surface, shaded grey->white, with black 500 m surface contours.
    ice_s = map((s, m) -> m ? s : NaN, surf, has_ice)
    heatmap!(ax, lon, lat, ice_s; colormap = ICE_CMAP, colorrange = ICE_RANGE,
             highclip = ICE_HIGH, interpolate = true)
    contour!(ax, lon, lat, ice_s; levels = 500:500:4500, color = :black, linewidth = 0.5)

    bsl = Float64(ds["bsl"][it])
    text!(ax, -176, -30; text = @sprintf("BSL = %.0f m", bsl),
          align = (:left, :center), fontsize = 13, color = :black)
    return nothing
end

function plot_run(label, path)
    ds = NCDataset(path)
    times = Float64.(ds["time"][:])
    lon, lat, pl, pt = grid_axes(ds)

    # panel index list: nearest slice to each requested time in TIMES_KA.
    its    = [nearest_index(times, t) for t in TIMES_KA]
    titles = [ (t == round(t) ? @sprintf("%.0f ka", t) : @sprintf("%.1f ka", t))
               for t in TIMES_KA ]

    fig = Figure(size = (1180, 760))
    Label(fig[0, 1:2], label; fontsize = 18, font = :bold)
    gmaps = fig[1, 1] = GridLayout()
    gbars = fig[1, 2] = GridLayout()

    for (j, it) in enumerate(its)
        r, c = fldmod1(j, 2)              # 2x2 layout: row, col
        ax = Axis(gmaps[r, c]; aspect = DataAspect(), title = titles[j],
                  xlabel = "", ylabel = "")
        limits!(ax, -180, 180, -90, 90)
        panel!(ax, ds, it, pl, pt, lon, lat)
    end

    # Three colorbars stacked on the right, each shrunk to ~2/3 of its cell.
    Colorbar(gbars[1, 1]; colormap = TOPO_CMAP, limits = TOPO_RANGE,
             lowclip = TOPO_LOW, highclip = TOPO_HIGH,
             label = "bedrock elevation [m]", height = Relative(0.66), valign = :center)
    Colorbar(gbars[2, 1]; colormap = ICE_CMAP, limits = ICE_RANGE,
             highclip = ICE_HIGH,
             label = "ice surface elevation [m]", height = Relative(0.66), valign = :center)
    Colorbar(gbars[3, 1]; colormap = RSL_CMAP, limits = RSL_CRANGE,
             ticks = RSL_CBTICKS, ticklabelsize = 9, lowclip = RSL_LOW, highclip = RSL_HIGH,
             label = "relative sea level [m]", height = Relative(0.95), valign = :center)
    rowsize!(gbars, 3, Auto(2.4))        # RSL bar needs room for its many labels
    colsize!(fig.layout, 2, Auto(0.18))
    close(ds)
    fn = joinpath(OUTDIR, "maps_" * replace(label, r"[^A-Za-z0-9]" => "_") * ".png")
    save(fn, fig); println("wrote ", fn)
end

function plot_bsl(runs)
    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time [ka]", ylabel = "barystatic sea level [m]",
              title = "Barystatic sea level vs present-day reference")
    for (label, path) in runs
        isfile(path) || (@warn "missing output, skipping BSL line" label path; continue)
        ds = NCDataset(path)
        t = Float64.(ds["time"][:]) ./ 1000
        b = Float64.(ds["bsl"][:])
        keep = t .>= BSL_TMIN_KA
        lines!(ax, t[keep], b[keep]; label = label, linewidth = 2)
        close(ds)
    end
    xlims!(ax, BSL_TMIN_KA, 0.0)
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
    plot_bsl(BSL_RUNS)
end

main()
