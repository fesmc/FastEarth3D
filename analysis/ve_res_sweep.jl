#!/usr/bin/env julia
# ve_res_sweep.jl — COMPUTE step for the VE-solver resolution + sub-step sweep
# staged by scripts/run_ve_res_sweep.sh. Quantifies how the FULL viscoelastic
# solver's answer and cost change with spherical-harmonic resolution (lmax) and the
# explicit sub-step ceiling (cfl), on the real 3-D deglaciation (Bagge 2021 lateral
# viscosity, Tarasov forcing). The production target lmax=128 / cfl=1 is the error
# reference; every other run is compared against it.
#
# Cross-resolution: runs at different lmax live on different Gauss grids, so each
# candidate's rsl is bilinearly regridded onto the reference (lmax=128) grid before
# the area-weighted error is formed (bsl is a global scalar — compared directly).
#
# Reads the raw run outputs and distils into a single portable file:
#   - analysis/ve_res_sweep_results.jld2   (grid, metrics, cost breakdown, bsl, resid maps)
#   - analysis/ve_res_sweep_summary.md     (human-readable table)
#
# Usage:  julia analysis/ve_res_sweep.jl [EXP_ROOT]
#         EXP_ROOT defaults to runs/ve_res_sweep

using NCDatasets
using JLD2
using Printf
using Statistics

const EXP       = length(ARGS) >= 1 ? ARGS[1] : "runs/ve_res_sweep"
const LMAX_LIST = [32, 64, 96, 128]
const LMAX_REF  = 128
const CFL_PROBE  = ["0.5", "1.5", "2.0", "2.5"] # extra cfl points at LMAX_REF (cfl=1 is the ref)
const VTOL_PROBE = ["0.1", "0.3", "1.0"]        # extra visc3d_tol [dex] at LMAX_REF (1e-3 is the ref)
const VTOL_REF   = 1.0e-3                        # reference 3-D-split threshold [dex]
const TIMES_KA  = [0.0, -10.0, -20.0, -26.0]    # residual-map snapshots [ka]
const JLD2_OUT  = "analysis/ve_res_sweep_results.jld2"
const MD_OUT    = "analysis/ve_res_sweep_summary.md"

# ---------------------------------------------------------------------------
# Candidate enumeration (mirrors run_ve_res_sweep.sh dir layout)
# ---------------------------------------------------------------------------
struct Cand
    label::String
    lmax::Int
    cfl::Float64
    vtol::Float64
    dir::String
end

function candidates()
    cs = Cand[]
    for L in LMAX_LIST
        push!(cs, Cand("lmax.$L", L, 1.0, VTOL_REF, joinpath(EXP, "lmax.$L")))
    end
    for c in CFL_PROBE
        push!(cs, Cand("lmax.$LMAX_REF.cfl$c", LMAX_REF, parse(Float64, c), VTOL_REF,
                       joinpath(EXP, "lmax.$LMAX_REF.cfl$c")))
    end
    for v in VTOL_PROBE
        push!(cs, Cand("lmax.$LMAX_REF.vtol$v", LMAX_REF, 1.0, parse(Float64, v),
                       joinpath(EXP, "lmax.$LMAX_REF.vtol$v")))
    end
    return cs
end

refdir() = joinpath(EXP, "lmax.$LMAX_REF")

# ---------------------------------------------------------------------------
# Cost: parse the [PROFILE] block of out.out (incl. the solid_earth_update split)
# ---------------------------------------------------------------------------
"Return Dict(se,drift,mem,read,write,nsolve,ne3d,ne) from out.out (NaN if absent)."
function read_profile(dir)
    f = joinpath(dir, "out.out")
    se = dr = mm = rd = wr = ns = n3 = ne = NaN
    if isfile(f)
        for ln in eachline(f)
            m = match(r"solid_earth_update\s*=\s*([\d.]+)\s*ms", ln); m !== nothing && (se = parse(Float64, m[1]))
            m = match(r"drift solve.*?=\s*([\d.]+)\s*ms", ln);        m !== nothing && (dr = parse(Float64, m[1]))
            m = match(r"memory advance\s*=\s*([\d.]+)\s*ms", ln);     m !== nothing && (mm = parse(Float64, m[1]))
            m = match(r"read_ice.*?=\s*([\d.]+)\s*ms", ln);           m !== nothing && (rd = parse(Float64, m[1]))
            m = match(r"fe_write_step.*?=\s*([\d.]+)\s*ms", ln);      m !== nothing && (wr = parse(Float64, m[1]))
            m = match(r"n_solve=\s*([\d.]+)", ln);                    m !== nothing && (ns = parse(Float64, m[1]))
            m = match(r"visc3d split:\s*(\d+)\s+of\s+(\d+)", ln)
            m !== nothing && (n3 = parse(Float64, m[1]); ne = parse(Float64, m[2]))
        end
    end
    return Dict("se" => se, "drift" => dr, "mem" => mm, "read" => rd, "write" => wr,
                "nsolve" => ns, "ne3d" => n3, "ne" => ne)
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

area_weights(lat, nlon) = repeat(reshape(cos.(deg2rad.(lat)), 1, :), nlon, 1)
wrmse(d, w) = sqrt(sum(w .* d .^ 2) / sum(w))
nearest_index(times_yr, t_ka) = argmin(abs.(times_yr ./ 1000 .- t_ka))

"""Bilinear regrid of fc (nlon_c, nlat_c) on (lon_c, lat_c) onto target (lon_r, lat_r).
lon is treated as a uniform periodic axis; lat is bracketed (clamped at the poles).
When the grids coincide this is the identity to round-off."""
function regrid(fc, lon_c, lat_c, lon_r, lat_r)
    nlc = length(lon_c)
    spx = sortperm(mod.(lon_c, 360.0)); lonc = mod.(lon_c, 360.0)[spx]
    spy = sortperm(lat_c);              latc = lat_c[spy]
    f = fc[spx, spy]                                   # ascending lon, ascending lat
    dlon = 360.0 / nlc
    fr = Array{Float64}(undef, length(lon_r), length(lat_r))
    for (jr, y) in enumerate(lat_r)
        yy = clamp(y, latc[1], latc[end])
        j2 = clamp(searchsortedfirst(latc, yy), 2, length(latc)); j1 = j2 - 1
        ty = (yy - latc[j1]) / (latc[j2] - latc[j1])
        for (ir, x) in enumerate(lon_r)
            t  = (mod(x, 360.0) - lonc[1]) / dlon
            i0 = floor(Int, t); tx = t - i0
            i1 = mod(i0,   nlc) + 1                     # periodic wrap
            i2 = mod(i0 + 1, nlc) + 1
            fr[ir, jr] = (1 - tx) * ((1 - ty) * f[i1, j1] + ty * f[i1, j2]) +
                         tx       * ((1 - ty) * f[i2, j1] + ty * f[i2, j2])
        end
    end
    return fr
end

# ---------------------------------------------------------------------------
# Per-candidate comparison against the reference (regridded), error + resid maps
# ---------------------------------------------------------------------------
"""Compare candidate against the reference. Candidate rsl is regridded onto the
reference grid (lon_r,lat_r raw order); errors are area-weighted (w, raw order);
resid[:,:,k] = rsl(cand)-rsl(ref) at TIMES_KA[k], plot-ready (pl,pt)."""
function compare(refnc, candnc, lon_r, lat_r, w, pl, pt, tidx)
    ref = NCDataset(refnc); cand = NCDataset(candnc)
    try
        nt   = length(ref["time"][:])
        lon_c = Float64.(cand["lon"][:]); lat_c = Float64.(cand["lat"][:])
        rr = ref["rsl"]; rc = cand["rsl"]
        resid = Array{Float32}(undef, length(pl), length(pt), length(tidx))

        sse = 0.0; wsum = 0.0; maxabs = 0.0; rsl_pd = NaN
        for it in 1:nt
            cg = regrid(Float64.(rc[:, :, it]), lon_c, lat_c, lon_r, lat_r)
            d  = cg .- Float64.(rr[:, :, it])
            sse += sum(w .* d .^ 2); wsum += sum(w)
            maxabs = max(maxabs, maximum(abs, d))
            it == nt && (rsl_pd = wrmse(d, w))
            k = findfirst(==(it), tidx)
            k !== nothing && (resid[:, :, k] = Float32.(d[pl, pt]))
        end

        br = Float64.(ref["bsl"][:]); bc = Float64.(cand["bsl"][:])
        err = Dict("rsl_rmse"    => sqrt(sse / wsum),
                   "rsl_rmse_pd" => rsl_pd,
                   "rsl_maxabs"  => maxabs,
                   "bsl_rmse"    => sqrt(mean((bc .- br) .^ 2)),
                   "bsl_pd_err"  => abs(bc[end] - br[end]))
        return err, bc, resid
    finally
        close(ref); close(cand)
    end
end

# ---------------------------------------------------------------------------
# Gather everything into a plain-Dict bundle (JLD2-friendly)
# ---------------------------------------------------------------------------
function gather()
    refnc = joinpath(refdir(), "out.nc")
    ref   = NCDataset(refnc)
    lon, lat, pl, pt = grid_axes(ref)
    lon_r = Float64.(ref["lon"][:]); lat_r = Float64.(ref["lat"][:])   # raw ref grid (regrid target)
    time  = Float64.(ref["time"][:])
    w     = area_weights(lat_r, length(lon_r))                        # raw-order weights
    tidx  = [nearest_index(time, t) for t in TIMES_KA]

    ref_bsl = Float64.(ref["bsl"][:])
    ref_rsl = Array{Float32}(undef, length(pl), length(pt), length(tidx))
    for (k, it) in enumerate(tidx)
        ref_rsl[:, :, k] = Float32.(Float64.(ref["rsl"][:, :, it])[pl, pt])
    end
    scales = Dict("rsl_pd_absmax" => maximum(abs, Float64.(ref["rsl"][:, :, end])),
                  "bsl_min" => minimum(ref_bsl), "bsl_pd" => ref_bsl[end])
    close(ref)
    ref_prof = read_profile(refdir())

    # Process each candidate, tolerating runs that are missing, still in progress
    # (fewer time slices than the reference), or otherwise unreadable — skip with a
    # warning so the analysis still distils whatever has finished.
    cands = Dict[];  skipped = String[];  nref = length(time)
    for c in candidates()
        nc = joinpath(c.dir, "out.nc")
        if !isfile(nc)
            push!(skipped, c.label);  @warn "no out.nc — skipping" label=c.label;  continue
        end
        try
            ntc = NCDataset(ds -> length(ds["time"][:]), nc)   # current slice count
            if ntc < nref
                push!(skipped, c.label)
                @warn "incomplete run — skipping" label=c.label slices="$ntc/$nref"
                continue
            end
            err, bsl, resid = compare(refnc, nc, lon_r, lat_r, w, pl, pt, tidx)
            prof = read_profile(c.dir)
            push!(cands, Dict("label" => c.label, "lmax" => c.lmax, "cfl" => c.cfl,
                              "vtol" => c.vtol, "err" => err, "prof" => prof,
                              "speedup" => ref_prof["se"] / prof["se"],
                              "bsl" => bsl, "resid" => resid))
            @printf("  %-18s rsl_rmse=%7.3f m  max=%7.2f m  cost=%7.1f ms  speedup=%5.2fx\n",
                    c.label, err["rsl_rmse"], err["rsl_maxabs"], prof["se"], ref_prof["se"] / prof["se"])
        catch e
            push!(skipped, c.label);  @warn "could not process — skipping" label=c.label exception=e
        end
    end
    isempty(skipped) || @info "skipped $(length(skipped)) run(s) (incomplete/missing): " * join(skipped, ", ")

    return Dict("lon" => lon, "lat" => lat, "time" => time, "times_ka" => TIMES_KA,
                "tidx" => tidx, "ref_prof" => ref_prof, "scales" => scales,
                "ref_bsl" => ref_bsl, "ref_rsl" => ref_rsl, "cands" => cands)
end

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
function table_md(S)
    io = IOBuffer()
    vp = S["ref_prof"]; sc = S["scales"]
    println(io, "# VE resolution & sub-step sweep — accuracy & cost\n")
    @printf(io, "Reference: lmax=%d, cfl=1. Every run is regridded onto the reference grid; ", LMAX_REF)
    println(io, "`rsl RMSE` is area-weighted (cos lat) over all space-time, `rsl PD` at present day.\n")
    @printf(io, "Reference cost: solid_earth_update = %.1f ms/step (drift %.1f + memory %.1f), n_solve = %.1f\n\n",
            vp["se"], vp["drift"], vp["mem"], vp["nsolve"])
    @printf(io, "Reference scales: |rsl|max(PD) = %.1f m, BSL min = %.1f m, BSL(PD) = %.2f m\n\n",
            sc["rsl_pd_absmax"], sc["bsl_min"], sc["bsl_pd"])
    println(io, "| run | lmax | cfl | vtol | ne3d | rsl RMSE | rsl PD | rsl max | bsl PD | cost | mem | nsolve | speedup |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in S["cands"]
        e = r["err"]; p = r["prof"]
        @printf(io, "| %s | %d | %.2g | %.0e | %.0f | %.3f m | %.3f m | %.3f m | %.3f m | %.1f ms | %.1f | %.1f | %.2fx |\n",
            r["label"], r["lmax"], r["cfl"], r["vtol"], p["ne3d"], e["rsl_rmse"], e["rsl_rmse_pd"], e["rsl_maxabs"],
            e["bsl_pd_err"], p["se"], p["mem"], p["nsolve"], r["speedup"])
    end
    println(io)
    return String(take!(io))
end

# ---------------------------------------------------------------------------
function main()
    println("== ve_res_sweep ==  exp=", EXP, "  ref=lmax.$LMAX_REF")
    S = gather()
    results = Dict{String,Any}("exp" => EXP, "ref_lmax" => LMAX_REF,
                               "lmax_list" => LMAX_LIST, "cfl_probe" => CFL_PROBE,
                               "times_ka" => TIMES_KA)
    for (k, v) in S
        results[k] = v
    end

    jldsave(JLD2_OUT; results = results)
    write(MD_OUT, table_md(S))
    sz = round(filesize(JLD2_OUT) / 1e6, digits = 1)
    println("wrote ", JLD2_OUT, "  (", sz, " MB)")
    println("wrote ", MD_OUT)
end

main()
