#!/usr/bin/env julia
# modal_vs_ve.jl — COMPUTE step. Quantify how well the reduced MODAL response
# approximates the full viscoelastic (VE) solver, in both accuracy and cost, for
# the experiment set staged by scripts/run_modal_vs_ve.sh.
#
# This reads the raw (large, ~3.7 GB) run outputs and distils everything needed
# for plotting into a single portable file:
#   - analysis/modal_vs_ve_results.jld2   (grid, metrics, bsl series, residual maps)
#   - analysis/modal_vs_ve_summary.md     (human-readable summary tables)
# Plot from the JLD2 with analysis/plot_modal_vs_ve.jl (no raw data needed).
#
# Two sets, each measured against the VE run in the *same* set (the ground truth):
#   radial   — 1-D (radial) viscosity. modal(n_modes=all) -> VE in the clean limit;
#              isolates the accuracy of the mode-count / ranking dial.
#   deglac3d — 3-D (laterally varying) viscosity, 1-D spin-up. Here modal is a
#              genuine approximation to VE.
#
# Usage:  julia analysis/modal_vs_ve.jl [EXP_ROOT]
#         EXP_ROOT defaults to runs/modal_vs_ve

using NCDatasets
using JLD2
using Printf
using Statistics

const EXP    = length(ARGS) >= 1 ? ARGS[1] : "runs/modal_vs_ve"
const SETS   = ["radial", "deglac3d"]
const NMODES = [1, 2, 4, 8]
const RANKS  = ["isostatic", "rate", "residue"]
const TIMES_KA = [0.0, -10.0, -20.0, -26.0]    # residual-map snapshots [ka]: PD, 10, 20, 26
const JLD2_OUT = "analysis/modal_vs_ve_results.jld2"
const MD_OUT   = "analysis/modal_vs_ve_summary.md"

# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------
struct Cand
    label::String   # short label
    nmodes::Int     # -1 == "all"
    rank::String    # "all" for modal_all
    dir::String
end

function candidates(set)
    root = joinpath(EXP, set)
    cs = Cand[Cand("all", -1, "all", joinpath(root, "modal_all"))]
    for n in NMODES, r in RANKS
        push!(cs, Cand("n$n/$r", n, r,
                       joinpath(root, "modal", "nmds.$n.mdrnk.$r")))
    end
    return cs
end

# ---------------------------------------------------------------------------
# Cost: parse the [PROFILE] block of out.out
# ---------------------------------------------------------------------------
"Return Dict(se,read,write,nsolve) ms/step (NaN if absent) from a run dir's out.out."
function read_profile(dir)
    f = joinpath(dir, "out.out")
    se = rd = wr = ns = NaN
    if isfile(f)
        for ln in eachline(f)
            m = match(r"solid_earth_update\s*=\s*([\d.]+)\s*ms", ln); m !== nothing && (se = parse(Float64, m[1]))
            m = match(r"read_ice.*?=\s*([\d.]+)\s*ms", ln);      m !== nothing && (rd = parse(Float64, m[1]))
            m = match(r"fe_write_step.*?=\s*([\d.]+)\s*ms", ln); m !== nothing && (wr = parse(Float64, m[1]))
            m = match(r"n_solve=\s*([\d.]+)", ln);               m !== nothing && (ns = parse(Float64, m[1]))
        end
    end
    return Dict("se" => se, "read" => rd, "write" => wr, "nsolve" => ns)
end

# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------
"Plot-ready ascending axes + permutations: lon rolled to [-180,180), lat ascending."
function grid_axes(ds)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    lon = map(x -> x > 180 ? x - 360 : x, lon)
    pl = sortperm(lon); pt = sortperm(lat)
    return lon[pl], lat[pt], pl, pt
end

"cos(lat) area weights as a (nlon, nlat) matrix matching stored field layout."
area_weights(lat, nlon) = repeat(reshape(cos.(deg2rad.(lat)), 1, :), nlon, 1)

wrmse(d, w) = sqrt(sum(w .* d .^ 2) / sum(w))

nearest_index(times_yr, t_ka) = argmin(abs.(times_yr ./ 1000 .- t_ka))

# ---------------------------------------------------------------------------
# Per-candidate comparison: error metrics + residual maps at the snapshot times
# ---------------------------------------------------------------------------
"""Compare candidate out.nc against VE reference. Returns (err::Dict, bsl::Vector,
resid::Array{Float32,3}) where resid[:,:,k] = rsl(modal)-rsl(VE) at TIMES_KA[k],
plot-ready (ascending lon/lat)."""
function compare(refnc, candnc, w, pl, pt, tidx)
    ref = NCDataset(refnc); cand = NCDataset(candnc)
    try
        nt = length(ref["time"][:])
        rr = ref["rsl"]; rc = cand["rsl"]
        nlon = length(pl); nlat = length(pt)
        resid = Array{Float32}(undef, nlon, nlat, length(tidx))

        sse = 0.0; wsum = 0.0; maxabs = 0.0; rsl_pd = NaN
        for it in 1:nt
            d = Float64.(rc[:, :, it]) .- Float64.(rr[:, :, it])
            sse += sum(w .* d .^ 2); wsum += sum(w)
            maxabs = max(maxabs, maximum(abs, d))
            it == nt && (rsl_pd = wrmse(d, w))
            k = findfirst(==(it), tidx)
            k !== nothing && (resid[:, :, k] = Float32.(d[pl, pt]))
        end

        br = Float64.(ref["bsl"][:]); bc = Float64.(cand["bsl"][:])
        err = Dict("rsl_rmse"    => sqrt(sse / wsum),     # pooled over all space-time
                   "rsl_rmse_pd" => rsl_pd,               # present day (final step)
                   "rsl_maxabs"  => maxabs,
                   "bsl_rmse"    => sqrt(mean((bc .- br) .^ 2)),
                   "bsl_pd_err"  => abs(bc[end] - br[end]))
        return err, bc, resid
    finally
        close(ref); close(cand)
    end
end

# ---------------------------------------------------------------------------
# Gather one set into a plain-Dict bundle (JLD2-friendly)
# ---------------------------------------------------------------------------
function gather(set)
    root   = joinpath(EXP, set)
    refnc  = joinpath(root, "ve", "out.nc")
    ref    = NCDataset(refnc)
    lon, lat, pl, pt = grid_axes(ref)
    time   = Float64.(ref["time"][:])
    w      = area_weights(Float64.(ref["lat"][:]), length(ref["lon"][:]))  # raw-order weights
    tidx   = [nearest_index(time, t) for t in TIMES_KA]

    ve_bsl  = Float64.(ref["bsl"][:])
    ve_rsl  = Array{Float32}(undef, length(pl), length(pt), length(tidx))
    for (k, it) in enumerate(tidx)
        ve_rsl[:, :, k] = Float32.(Float64.(ref["rsl"][:, :, it])[pl, pt])
    end
    scales = Dict("rsl_pd_absmax" => maximum(abs, Float64.(ref["rsl"][:, :, end])),
                  "bsl_min" => minimum(ve_bsl), "bsl_pd" => ve_bsl[end])
    close(ref)
    ve_prof = read_profile(joinpath(root, "ve"))

    cands = Dict[]
    for c in candidates(set)
        nc = joinpath(c.dir, "out.nc")
        isfile(nc) || (@warn "missing out.nc, skipping" dir=c.dir; continue)
        err, bsl, resid = compare(refnc, nc, w, pl, pt, tidx)
        prof = read_profile(c.dir)
        push!(cands, Dict("label" => c.label, "nmodes" => c.nmodes, "rank" => c.rank,
                          "err" => err, "prof" => prof, "speedup" => ve_prof["se"] / prof["se"],
                          "bsl" => bsl, "resid" => resid))
        @printf("  %-12s rsl_rmse=%7.3f m  cost=%6.1f ms  speedup=%6.1fx\n",
                c.label, err["rsl_rmse"], prof["se"], ve_prof["se"] / prof["se"])
    end

    return Dict("lon" => lon, "lat" => lat, "time" => time, "times_ka" => TIMES_KA,
                "tidx" => tidx, "ve_prof" => ve_prof, "scales" => scales,
                "ve_bsl" => ve_bsl, "ve_rsl" => ve_rsl, "cands" => cands)
end

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
function table_md(set, S)
    io = IOBuffer()
    vp = S["ve_prof"]; sc = S["scales"]
    println(io, "### $set\n")
    @printf(io, "VE reference: solver cost = %.1f ms/step, n_solve = %.1f sub-steps/step\n\n",
            vp["se"], vp["nsolve"])
    @printf(io, "VE scales: |rsl|max(PD) = %.1f m, BSL min = %.1f m, BSL(PD) = %.2f m\n\n",
            sc["rsl_pd_absmax"], sc["bsl_min"], sc["bsl_pd"])
    println(io, "| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for r in S["cands"]
        e = r["err"]; p = r["prof"]
        @printf(io, "| %s | %.3f m | %.3f m | %.3f m | %.3f m | %.3f m | %.1f ms | %.1fx | %.1f |\n",
            r["label"], e["rsl_rmse"], e["rsl_rmse_pd"], e["rsl_maxabs"],
            e["bsl_rmse"], e["bsl_pd_err"], p["se"], r["speedup"], p["nsolve"])
    end
    println(io)
    return String(take!(io))
end

# ---------------------------------------------------------------------------
function main()
    results = Dict{String,Any}("sets" => SETS, "nmodes" => NMODES, "ranks" => RANKS,
                               "times_ka" => TIMES_KA, "exp" => EXP)
    md = IOBuffer()
    println(md, "# Modal vs VE — accuracy & cost summary\n")
    println(md, "Each modal run is compared against the VE run in the same set. ",
                "`rsl RMSE` is area-weighted (cos lat) over all space-time; `rsl PD` at ",
                "present day; `cost` is the per-coupling-step solver time; `speedup` = VE/modal.\n")
    for set in SETS
        println("== $set ==")
        S = gather(set)
        results[set] = S
        print(md, table_md(set, S))
    end

    jldsave(JLD2_OUT; results = results)
    write(MD_OUT, String(take!(md)))
    sz = round(filesize(JLD2_OUT) / 1e6, digits = 1)
    println("wrote ", JLD2_OUT, "  (", sz, " MB)")
    println("wrote ", MD_OUT)
end

main()
