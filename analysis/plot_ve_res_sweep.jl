#!/usr/bin/env julia
# plot_ve_res_sweep.jl — PLOT step for the VE-solver resolution + sub-step sweep.
# Portable: reads only the distilled JLD2 written by analysis/ve_res_sweep.jl
# (no raw run data needed), so it runs anywhere.
#
# Every run is compared against the lmax=128 / cfl=1 reference (regridded onto its
# grid by the compute step). Figures (analysis/figs/ve_res_sweep/):
#   accuracy_cost_vs_lmax.png  rsl error vs lmax + cost(total/memory/drift) vs lmax
#   pareto.png                 accuracy-cost trade-off (rsl RMSE vs solver ms/step)
#   cost_breakdown.png         drift / memory advance / SLE+rest, ms/step, per run
#   bsl_timeseries.png         reference vs each run, barystatic sea level
#   resid_maps_<t>ka.png       run-minus-reference rsl residual maps, all runs
#   rsl_panels_<t>ka.png       rows ref/l64/l32: col1 absolute rsl, col2 anomaly vs ref
#
# Usage:  julia analysis/plot_ve_res_sweep.jl [results.jld2]
#         defaults to analysis/ve_res_sweep_results.jld2

using JLD2
using CairoMakie
using GeoMakie   # thin continental outlines (110m Natural Earth, bundled offline)
using Printf

const JLD2_IN = length(ARGS) >= 1 ? ARGS[1] : "analysis/ve_res_sweep_results.jld2"
const OUTDIR  = "analysis/figs/ve_res_sweep"

# Distinct colour per run label (resolution sweep + cfl + visc3d_tol probes).
const PALETTE = [:firebrick, :darkorange, :goldenrod, :olive, :seagreen, :teal,
                 :royalblue, :navy, :purple, :magenta, :brown, :black]
runcolor(i) = PALETTE[mod1(i, length(PALETTE))]
const VTOL_REF = 1e-3                            # default 3-D-split threshold [dex]
vtol(c) = get(c, "vtol", VTOL_REF)               # tolerate older jld2 without the key

cand_by_label(C, lbl) = (i = findfirst(c -> c["label"] == lbl, C); i === nothing ? nothing : C[i])
# resolution ladder: cfl=1 AND default vtol (so the cfl/vtol probes don't leak in)
res_sweep(C)  = sort(filter(c -> c["cfl"] == 1.0 && vtol(c) == VTOL_REF, C), by = c -> c["lmax"])
cfl_probes(C, ref)  = sort(filter(c -> c["lmax"] == ref && c["cfl"] != 1.0, C), by = c -> c["cfl"])
vtol_probes(C, ref) = sort(filter(c -> c["lmax"] == ref && c["cfl"] == 1.0 && vtol(c) != VTOL_REF, C),
                           by = c -> vtol(c))
rest(p) = max(p["se"] - p["drift"] - p["mem"], 0.0)                              # SLE + coupling
efloor(x) = max(x, 1e-4)                                                         # log-axis floor

# ---------------------------------------------------------------------------
# fig 1: accuracy and cost vs lmax (the resolution ladder, cfl=1).
function fig_accuracy_cost(R)
    C = R["cands"]; sweep = res_sweep(C)
    ls   = [Float64(c["lmax"])      for c in sweep]
    cost = [c["prof"]["se"]/1000    for c in sweep]   # ms -> s
    mem  = [c["prof"]["mem"]/1000   for c in sweep]
    dr   = [c["prof"]["drift"]/1000 for c in sweep]
    # accuracy: drop the reference itself (lmax_ref) — its self-error is exactly 0.
    acc  = filter(c -> c["lmax"] != R["ref_lmax"], sweep)
    als  = [Float64(c["lmax"])             for c in acc]
    rmse = [efloor(c["err"]["rsl_rmse"])   for c in acc]
    mx   = [efloor(c["err"]["rsl_maxabs"]) for c in acc]

    fig = Figure(size = (1100, 470))
    ax1 = Axis(fig[1, 1]; xlabel = "lmax", ylabel = "rsl error vs lmax$(R["ref_lmax"]) [m]",
               title = "accuracy (vs lmax$(R["ref_lmax"]) reference)", xscale = log2, yscale = log10,
               xticks = (ls, string.(Int.(ls))))
    scatterlines!(ax1, als, rmse; color = :seagreen, marker = :circle, markersize = 11,
                  linewidth = 2, label = "area-weighted RMSE")
    scatterlines!(ax1, als, mx; color = :firebrick, marker = :rect, markersize = 11,
                  linewidth = 2, linestyle = :dash, label = "max |Δrsl|")
    axislegend(ax1; position = :rt, framevisible = false)

    ax2 = Axis(fig[1, 2]; xlabel = "lmax", ylabel = "cost [s/step]",
               title = "cost", xscale = log2,
               xticks = (ls, string.(Int.(ls))))
    scatterlines!(ax2, ls, cost; color = :black,    marker = :circle, markersize = 11,
                  linewidth = 2, label = "solid_earth_update")
    scatterlines!(ax2, ls, mem;  color = :royalblue, marker = :utriangle, markersize = 11,
                  linewidth = 2, label = "memory advance")
    scatterlines!(ax2, ls, dr;   color = :darkorange, marker = :dtriangle, markersize = 11,
                  linewidth = 2, label = "drift solve")
    axislegend(ax2; position = :lt, framevisible = false)

    Label(fig[0, 1:2], "VE solver vs resolution (reference = lmax$(R["ref_lmax"]), cfl=1)";
          fontsize = 16, font = :bold)
    fn = joinpath(OUTDIR, "accuracy_cost_vs_lmax.png"); save(fn, fig); println("wrote ", fn)
end

# fig 2: Pareto — cost (ms/step) vs rsl RMSE, every run. Reference cost marked.
function fig_pareto(R)
    C = R["cands"]
    fig = Figure(size = (760, 560))
    ax = Axis(fig[1, 1]; xlabel = "solver cost [s/step]", ylabel = "rsl RMSE vs lmax$(R["ref_lmax"]) [m]",
              title = "accuracy–cost trade-off (down-left is better)",
              yscale = log10, xticks = 0:3:15)   # linear cost axis (plain numbers), log RMSE
    refcost = R["ref_prof"]["se"]/1000   # ms -> s
    vlines!(ax, [refcost]; color = :gray, linestyle = :dash)
    text!(ax, refcost, 1.0; text = @sprintf("ref = %.1f s", refcost),   # up into empty space, nudged left
          align = (:center, :center), space = :data, offset = (-9, 0), rotation = pi/2,
          fontsize = 10, color = :gray40)
    # dots in candidate order (stable colours)
    for (i, c) in enumerate(C)
        scatter!(ax, [c["prof"]["se"]/1000], [efloor(c["err"]["rsl_rmse"])];
                 color = runcolor(i), markersize = 14)
    end
    # labels in cost order with an alternating vertical nudge so points at a similar RMSE
    # (e.g. cfl1.5 / cfl2.0 / vtol1.0 ~0.08 m) don't overprint; drop the common "lmax." prefix.
    for (rank, i) in enumerate(sortperm([c["prof"]["se"] for c in C]))
        c = C[i]
        text!(ax, c["prof"]["se"]/1000, efloor(c["err"]["rsl_rmse"]);
              text = replace(c["label"], "lmax." => ""), align = (:left, :center),
              offset = (8, isodd(rank) ? 9 : -9), fontsize = 10, color = runcolor(i))
    end
    # pad both axes so the right-side / edge labels fit
    costs = [c["prof"]["se"]/1000 for c in C]
    errs  = [efloor(c["err"]["rsl_rmse"]) for c in C]
    xlims!(ax, 0, 15)                                # matches xticks 0:3:15
    ylims!(ax, minimum(errs) * 0.4,  maximum(errs) * 3)
    fn = joinpath(OUTDIR, "pareto.png"); save(fn, fig); println("wrote ", fn)
end

# fig 3: cost breakdown — drift / memory advance / SLE+rest, ms/step, per run.
function fig_cost_breakdown(R)
    C = R["cands"]
    labels = [c["label"] for c in C]
    x = 1:length(C)
    dr   = [c["prof"]["drift"]/1000 for c in C]   # ms -> s
    mem  = [c["prof"]["mem"]/1000   for c in C]
    rst  = [rest(c["prof"])/1000    for c in C]
    # stacked bars: drift (bottom) + memory + rest
    fig = Figure(size = (260 + 120 * length(C), 520))
    ax = Axis(fig[1, 1]; ylabel = "s / coupling step", title = "solid_earth_update cost breakdown",
              xticks = (collect(x), labels), xticklabelrotation = pi/6)
    cols = [:darkorange, :royalblue, :gray70]
    barplot!(ax, repeat(collect(x), 3),
             vcat(dr, mem, rst);
             stack = vcat(fill(1, length(C)), fill(2, length(C)), fill(3, length(C))),
             color = vcat(fill(cols[1], length(C)), fill(cols[2], length(C)), fill(cols[3], length(C))))
    elems = [PolyElement(color = c) for c in cols]
    Legend(fig[1, 2], elems, ["drift solve (band LU)", "memory advance (3-D SHT)", "SLE + coupling"];
           framevisible = false)
    fn = joinpath(OUTDIR, "cost_breakdown.png"); save(fn, fig); println("wrote ", fn)
end

# fig 4: BSL time series — reference vs every run.
function fig_bsl(R)
    C = R["cands"]; t = R["time"] ./ 1000
    fig = Figure(size = (900, 520))
    ax = Axis(fig[1, 1]; xlabel = "time [ka]", ylabel = "BSL [m]",
              title = "barystatic sea level: reference vs each run")
    lines!(ax, t[2:end], R["ref_bsl"][2:end]; color = :black, linewidth = 3, label = "ref lmax$(R["ref_lmax"])")
    for (i, c) in enumerate(C)
        c["label"] == "lmax.$(R["ref_lmax"])" && continue   # identical to ref
        lines!(ax, t[2:end], c["bsl"][2:end]; color = runcolor(i), linewidth = 1.8,
               linestyle = :dash, label = c["label"])
    end
    axislegend(ax; position = :rb, framevisible = false)
    fn = joinpath(OUTDIR, "bsl_timeseries.png"); save(fn, fig); println("wrote ", fn)
end

# fig 5: residual maps — one figure per snapshot time, all runs.
function fig_resid_maps(R)
    C = R["cands"]; lon = R["lon"]; lat = R["lat"]
    for (k, tka) in enumerate(R["times_ka"])
        # robust symmetric range: 99th pct of |resid| across runs at this time, ≥5 m.
        allabs = sort(vcat([vec(abs.(Float64.(c["resid"][:, :, k]))) for c in C]...))
        p99 = allabs[clamp(round(Int, 0.99 * length(allabs)), 1, length(allabs))]
        cap = max(5.0, ceil(p99 / 5) * 5); crange = (-cap, cap)
        ncol = 3; nrow = cld(length(C), ncol)
        fig = Figure(size = (380 * ncol + 120, 240 * nrow + 70))
        local hm
        for (i, c) in enumerate(C)
            d = Float64.(c["resid"][:, :, k])
            ax = Axis(fig[cld(i, ncol), mod1(i, ncol)]; aspect = DataAspect(),
                      title = @sprintf("%s  (max %.1f m)", c["label"], maximum(abs, d)),
                      titlesize = 11)
            limits!(ax, -180, 180, -90, 90); hidedecorations!(ax)
            hm = heatmap!(ax, lon, lat, d; colormap = :balance, colorrange = crange)
        end
        Colorbar(fig[1:nrow, ncol + 1], hm; label = "rsl(run) − rsl(ref) [m]")
        Label(fig[0, 1:ncol], @sprintf("residual maps vs lmax%d at %.0f ka", R["ref_lmax"], tka);
              fontsize = 15, font = :bold)
        tag = @sprintf("%02dka", round(Int, abs(tka)))
        fn = joinpath(OUTDIR, "resid_maps_$(tag).png"); save(fn, fig); println("wrote ", fn)
    end
end

# fig 6: rsl panels — one figure per snapshot. Rows ref / l64 / l32 (cfl=1 ladder).
# Col 1 = absolute rsl, col 2 = anomaly vs ref. Absolute = ref_rsl + resid.
function fig_rsl_panels(R)
    C = R["cands"]; lon = R["lon"]; lat = R["lat"]; ref = R["ref_lmax"]
    sweep = reverse(res_sweep(C))                       # ref (highest lmax) first
    rows = [(c["label"] == "lmax.$ref" ? "ref l$ref" : "l$(c["lmax"])",
             c["label"] == "lmax.$ref" ? nothing : c) for c in sweep]
    anom_range = (-30.0, 30.0)
    nrow = length(rows)
    coast = GeoMakie.coastlines()
    function style!(ax; bottom = false, xlabel = "")
        limits!(ax, -180, 180, -90, 90)
        lines!(ax, coast; color = (:gray20, 0.6), linewidth = 0.35)
        if bottom
            ax.xlabel = xlabel; ax.xlabelsize = 15
            hidexdecorations!(ax; label = false); hideydecorations!(ax)
        else
            hidedecorations!(ax)
        end
    end
    for (k, tka) in enumerate(R["times_ka"])
        refrsl = Float64.(R["ref_rsl"][:, :, k])
        vs = sort(vec(abs.(refrsl)))
        p99 = vs[clamp(round(Int, 0.99 * length(vs)), 1, length(vs))]
        rmax = max(50.0, ceil(p99 / 50) * 50); rsl_range = (-rmax, rmax)

        fig = Figure(size = (1000, 250 * nrow + 120))
        local hm_rsl, hm_anom
        for (i, (lbl, c)) in enumerate(rows)
            Label(fig[i, 0], lbl; font = :bold, fontsize = 16, rotation = pi/2)
            bottom = i == nrow
            af = c === nothing ? refrsl : refrsl .+ Float64.(c["resid"][:, :, k])
            ax1 = Axis(fig[i, 1]; titlesize = 14, title = @sprintf("max %.0f m", maximum(abs, af)))
            hm_rsl = heatmap!(ax1, lon, lat, af; colormap = :vik, colorrange = rsl_range)
            style!(ax1; bottom = bottom, xlabel = "absolute rsl")
            c === nothing && continue
            d = Float64.(c["resid"][:, :, k])
            ax2 = Axis(fig[i, 2]; titlesize = 14, title = @sprintf("max %.1f m", maximum(abs, d)))
            hm_anom = heatmap!(ax2, lon, lat, d; colormap = :balance, colorrange = anom_range)
            style!(ax2; bottom = bottom, xlabel = "anomaly vs ref")
        end
        for i in 1:nrow
            rowsize!(fig.layout, i, Aspect(1, 0.5))
        end
        Colorbar(fig[2:nrow, 3], hm_rsl;  label = @sprintf("rsl [m]  (clipped at ±%.0f)", rmax),
                 width = 11, height = Relative(0.42), labelsize = 12, ticklabelsize = 11)
        Colorbar(fig[2:nrow, 4], hm_anom; label = "rsl(run) − rsl(ref) [m]",
                 width = 11, height = Relative(0.42), labelsize = 12, ticklabelsize = 11)
        Label(fig[0, 0:4], @sprintf("rsl at %.0f ka  (resolution ladder vs lmax%d)", tka, ref);
              fontsize = 19, font = :bold)
        resize_to_layout!(fig)
        tag = @sprintf("%02dka", round(Int, abs(tka)))
        fn = joinpath(OUTDIR, "rsl_panels_$(tag).png"); save(fn, fig); println("wrote ", fn)
    end
end

# fig 7: visc3d_tol axis — accuracy + cost (with ne3d) vs the 3-D-split threshold.
function fig_visc3d_tol(R)
    C = R["cands"]; ref = R["ref_lmax"]
    probes = vtol_probes(C, ref)
    isempty(probes) && return                         # no vtol axis (e.g. older jld2)
    refc = cand_by_label(C, "lmax.$ref")
    sweep = sort(refc === nothing ? probes : vcat([refc], probes), by = c -> vtol(c))
    vs   = [vtol(c)             for c in sweep]
    cost = [c["prof"]["se"]/1000  for c in sweep]
    mem  = [c["prof"]["mem"]/1000 for c in sweep]
    ne3d = [c["prof"]["ne3d"]     for c in sweep]
    acc  = filter(c -> vtol(c) != VTOL_REF, sweep)    # drop the reference (error 0)
    avs  = [vtol(c)                        for c in acc]
    rmse = [efloor(c["err"]["rsl_rmse"])   for c in acc]
    mx   = [efloor(c["err"]["rsl_maxabs"]) for c in acc]

    fig = Figure(size = (1100, 470))
    ax1 = Axis(fig[1, 1]; xlabel = "visc3d_tol [dex]", ylabel = "rsl error vs reference [m]",
               title = "accuracy", xscale = log10, yscale = log10)
    scatterlines!(ax1, avs, rmse; color = :seagreen, marker = :circle, markersize = 11,
                  linewidth = 2, label = "area-weighted RMSE")
    scatterlines!(ax1, avs, mx; color = :firebrick, marker = :rect, markersize = 11,
                  linewidth = 2, linestyle = :dash, label = "max |Δrsl|")
    axislegend(ax1; position = :lt, framevisible = false)

    ax2 = Axis(fig[1, 2]; xlabel = "visc3d_tol [dex]", ylabel = "cost [s/step]",
               title = "cost (ne3d labelled)", xscale = log10)
    scatterlines!(ax2, vs, cost; color = :black, marker = :circle, markersize = 11,
                  linewidth = 2, label = "solid_earth_update")
    scatterlines!(ax2, vs, mem; color = :royalblue, marker = :utriangle, markersize = 11,
                  linewidth = 2, label = "memory advance")
    for (x, y, n) in zip(vs, cost, ne3d)
        text!(ax2, x, y; text = @sprintf("ne3d=%.0f", n), offset = (6, 4), fontsize = 9, color = :gray30)
    end
    axislegend(ax2; position = :rt, framevisible = false)

    Label(fig[0, 1:2], "VE solver vs visc3d_tol (lmax$ref, cfl=1; higher tol → fewer 3-D elements)";
          fontsize = 16, font = :bold)
    fn = joinpath(OUTDIR, "visc3d_tol.png"); save(fn, fig); println("wrote ", fn)
end

# Drop runs that went numerically unstable (e.g. cfl past the explicit M≤2 limit):
# their rsl/bsl carry absurd values that would wreck every shared axis. Judged against
# the physical scale of the reference solution; the summary table still records them.
function drop_unstable!(R)
    phys = 10 * R["scales"]["rsl_pd_absmax"]            # |value| beyond this ⇒ diverged
    ok(x) = isfinite(x) && abs(x) <= phys
    stable(c) = ok(c["err"]["rsl_rmse"]) && ok(c["err"]["rsl_maxabs"]) && all(ok, c["bsl"])
    bad = filter(c -> !stable(c), R["cands"])
    isempty(bad) || @info "excluding unstable run(s) from all figures: " *
                          join([c["label"] for c in bad], ", ")
    R["cands"] = filter(stable, R["cands"])
    return R
end

# ---------------------------------------------------------------------------
function main()
    mkpath(OUTDIR)
    R = load(JLD2_IN, "results")
    drop_unstable!(R)
    fig_accuracy_cost(R)
    fig_pareto(R)
    fig_cost_breakdown(R)
    fig_visc3d_tol(R)
    fig_bsl(R)
    fig_resid_maps(R)
    fig_rsl_panels(R)
end

main()
