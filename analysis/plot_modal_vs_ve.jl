#!/usr/bin/env julia
# plot_modal_vs_ve.jl — PLOT step. Portable: reads only the distilled JLD2 written
# by analysis/modal_vs_ve.jl (no raw run data needed), so it runs anywhere.
#
# Figures (analysis/figs/modal_vs_ve/):
#   convergence.png         rsl RMSE vs n_modes, a line per mode_rank, per set
#   pareto.png              accuracy-cost trade-off (rsl RMSE vs solver ms/step)
#   bsl_timeseries.png      VE vs cheapest modal (n=1) + all, barystatic sea level
#   resid_maps_<set>_<t>ka  modal-VE residual maps at PD/10/20/26 ka, all candidates
#   rsl_panels_<set>_<t>ka  one figure per time slice: rows VE/modal-all/8/4/1,
#                           col 1 = absolute rsl, col 2 = anomaly vs VE
#
# Usage:  julia analysis/plot_modal_vs_ve.jl [results.jld2] [rsl_panel_set]
#         defaults to analysis/modal_vs_ve_results.jld2  deglac3d

using JLD2
using CairoMakie
using GeoMakie   # for thin continental outlines (110m Natural Earth, bundled offline)
using Printf

const JLD2_IN  = length(ARGS) >= 1 ? ARGS[1] : "analysis/modal_vs_ve_results.jld2"
const RSL_SET  = length(ARGS) >= 2 ? ARGS[2] : "deglac3d"   # set for the rsl-panel figures
const RSL_RANK = "rate"                                     # mode ranking used for modal-8/4/1
const OUTDIR   = "analysis/figs/modal_vs_ve"
const RANK_COLOR = Dict("isostatic" => :firebrick, "rate" => :seagreen,
                        "residue" => :royalblue, "all" => :black)

cand_by_label(S, lbl) = (i = findfirst(c -> c["label"] == lbl, S["cands"]);
                         i === nothing ? nothing : S["cands"][i])
cands_by_rank(S, r) = filter(c -> c["rank"] == r, S["cands"])

# ---------------------------------------------------------------------------
# fig 1: convergence — rsl RMSE vs n_modes, a line per rank, one axis per set.
function fig_convergence(R, sets)
    fig = Figure(size = (1100, 460))
    for (j, set) in enumerate(sets)
        S = R[set]
        ax = Axis(fig[1, j]; xlabel = "n_modes per degree", ylabel = "rsl RMSE vs VE [m]",
                  title = set, yscale = log10, xscale = log2,
                  xticks = (R["nmodes"], string.(R["nmodes"])))
        for r in R["ranks"]
            cs = sort(cands_by_rank(S, r), by = c -> c["nmodes"])
            xs = [Float64(c["nmodes"]) for c in cs]
            ys = [max(c["err"]["rsl_rmse"], 1e-6) for c in cs]
            scatterlines!(ax, xs, ys; color = RANK_COLOR[r], label = r,
                          marker = :circle, markersize = 11, linewidth = 2)
        end
        allc = cand_by_label(S, "all")
        allc !== nothing && hlines!(ax, [max(allc["err"]["rsl_rmse"], 1e-6)];
                                    color = :black, linestyle = :dash, label = "all")
        j == 1 && axislegend(ax; position = :rt, framevisible = false)
    end
    Label(fig[0, 1:2], "Modal accuracy vs mode count (lower = closer to VE)";
          fontsize = 16, font = :bold)
    fn = joinpath(OUTDIR, "convergence.png"); save(fn, fig); println("wrote ", fn)
end

# fig 2: Pareto — cost (ms/step) vs rsl RMSE, one axis per set. VE cost as vline.
function fig_pareto(R, sets)
    fig = Figure(size = (1100, 480))
    for (j, set) in enumerate(sets)
        S = R[set]
        ax = Axis(fig[1, j]; xlabel = "solver cost [ms/step]", ylabel = "rsl RMSE vs VE [m]",
                  title = set, xscale = log10, yscale = log10)
        vlines!(ax, [S["ve_prof"]["se"]]; color = :gray, linestyle = :dash)
        text!(ax, S["ve_prof"]["se"], 0.0;
              text = @sprintf(" VE = %.0f ms", S["ve_prof"]["se"]),
              align = (:left, :bottom), space = :data, offset = (2, 2),
              rotation = pi/2, fontsize = 10, color = :gray40)
        for r in vcat(R["ranks"], "all")
            cs = cands_by_rank(S, r)
            isempty(cs) && continue
            xs = [c["prof"]["se"] for c in cs]
            ys = [max(c["err"]["rsl_rmse"], 1e-6) for c in cs]
            scatter!(ax, xs, ys; color = RANK_COLOR[r], markersize = 13, label = r)
            for c in cs
                text!(ax, c["prof"]["se"], max(c["err"]["rsl_rmse"], 1e-6);
                      text = c["nmodes"] == -1 ? "all" : string(c["nmodes"]),
                      align = (:left, :center), offset = (6, 0), fontsize = 9,
                      color = RANK_COLOR[r])
            end
        end
        j == 1 && axislegend(ax; position = :lb, framevisible = false)
    end
    Label(fig[0, 1:2], "Accuracy–cost trade-off (down-left is better; VE cost marked)";
          fontsize = 16, font = :bold)
    fn = joinpath(OUTDIR, "pareto.png"); save(fn, fig); println("wrote ", fn)
end

# fig 3: BSL time series — VE vs n_modes=all and the cheapest (n=1) per rank.
function fig_bsl(R, sets)
    fig = Figure(size = (1100, 460))
    for (j, set) in enumerate(sets)
        S = R[set]
        t = S["time"] ./ 1000
        ax = Axis(fig[1, j]; xlabel = "time [ka]", ylabel = "BSL [m]", title = set)
        lines!(ax, t[2:end], S["ve_bsl"][2:end]; color = :black, linewidth = 3, label = "VE")
        for lbl in ["all", "n1/isostatic", "n1/rate", "n1/residue"]
            c = cand_by_label(S, lbl); c === nothing && continue
            col = c["rank"] == "all" ? :black : RANK_COLOR[c["rank"]]
            lines!(ax, t[2:end], c["bsl"][2:end]; color = col,
                   linestyle = c["rank"] == "all" ? :dot : :dash,
                   linewidth = 2, label = lbl)
        end
        j == 2 && axislegend(ax; position = :rb, framevisible = false)
    end
    Label(fig[0, 1:2], "Barystatic sea level: VE vs cheapest modal (n_modes=1) + all";
          fontsize = 16, font = :bold)
    fn = joinpath(OUTDIR, "bsl_timeseries.png"); save(fn, fig); println("wrote ", fn)
end

# fig 4: residual maps — one figure per (set, snapshot time), all 13 candidates.
function fig_resid_maps(R, sets)
    crange = (-30.0, 30.0)
    for set in sets
        S = R[set]
        lon = S["lon"]; lat = S["lat"]
        for (k, tka) in enumerate(R["times_ka"])
            cs = S["cands"]
            ncand = length(cs); ncol = 4; nrow = cld(ncand, ncol)
            fig = Figure(size = (360 * ncol + 120, 230 * nrow + 70))
            local hm
            for (i, c) in enumerate(cs)
                d = Float64.(c["resid"][:, :, k])
                ax = Axis(fig[cld(i, ncol), mod1(i, ncol)]; aspect = DataAspect(),
                          title = @sprintf("%s  (max %.1f m)", c["label"], maximum(abs, d)),
                          titlesize = 11)
                limits!(ax, -180, 180, -90, 90); hidedecorations!(ax)
                hm = heatmap!(ax, lon, lat, d; colormap = :balance, colorrange = crange)
            end
            Colorbar(fig[1:nrow, ncol + 1], hm; label = "rsl(modal) − rsl(VE) [m]")
            Label(fig[0, 1:ncol], @sprintf("%s: residual maps at %.0f ka", set, tka);
                  fontsize = 15, font = :bold)
            tag = @sprintf("%02dka", round(Int, abs(tka)))
            fn = joinpath(OUTDIR, "resid_maps_$(set)_$(tag).png")
            save(fn, fig); println("wrote ", fn)
        end
    end
end

# fig 5: rsl panels — one figure per snapshot time. Rows VE / modal-all / 8 / 4 / 1
# (modal-N uses the `rank` ranking). Col 1 = absolute rsl, col 2 = anomaly vs VE.
# modal absolute is reconstructed as ve_rsl + resid (resid = rsl_modal − rsl_VE).
function fig_rsl_panels(R, set; rank = RSL_RANK)
    S = R[set]; lon = S["lon"]; lat = S["lat"]
    rows = [("VE",        nothing),
            ("modal-all", cand_by_label(S, "all")),
            ("modal-8",   cand_by_label(S, "n8/$rank")),
            ("modal-4",   cand_by_label(S, "n4/$rank")),
            ("modal-1",   cand_by_label(S, "n1/$rank"))]
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
        ve = Float64.(S["ve_rsl"][:, :, k])
        absflds = [c === nothing ? ve : ve .+ Float64.(c["resid"][:, :, k]) for (_, c) in rows]
        # shared rsl range: robust 99th-percentile of |VE| (raw max is one localised
        # ice-margin extreme that would wash out all structure), rounded up to 50 m.
        vs = sort(vec(abs.(ve)))
        p99 = vs[clamp(round(Int, 0.99 * length(vs)), 1, length(vs))]
        rmax = max(50.0, ceil(p99 / 50) * 50)
        rsl_range = (-rmax, rmax)

        fig = Figure(size = (1000, 250 * nrow + 120))
        local hm_rsl, hm_anom
        for (i, (lbl, c)) in enumerate(rows)
            Label(fig[i, 0], lbl; font = :bold, fontsize = 16, rotation = pi/2)
            bottom = i == nrow
            af = absflds[i]
            ax1 = Axis(fig[i, 1]; titlesize = 14,
                       title = @sprintf("max %.0f m", maximum(abs, af)))
            hm_rsl = heatmap!(ax1, lon, lat, af; colormap = :vik, colorrange = rsl_range)
            style!(ax1; bottom = bottom, xlabel = "absolute rsl")

            c === nothing && continue
            d = Float64.(c["resid"][:, :, k])
            ax2 = Axis(fig[i, 2]; titlesize = 14,
                       title = @sprintf("max %.1f m", maximum(abs, d)))
            hm_anom = heatmap!(ax2, lon, lat, d; colormap = :balance, colorrange = anom_range)
            style!(ax2; bottom = bottom, xlabel = "anomaly vs VE")
        end
        # tie every content row's height to the map-column width (lon:lat = 2:1) so
        # all panels — including the VE row that lacks an anomaly axis — match.
        for i in 1:nrow
            rowsize!(fig.layout, i, Aspect(1, 0.5))
        end
        Colorbar(fig[2:nrow, 3], hm_rsl;  label = @sprintf("rsl [m]  (clipped at ±%.0f)", rmax),
                 width = 11, height = Relative(0.42), labelsize = 12, ticklabelsize = 11)
        Colorbar(fig[2:nrow, 4], hm_anom; label = "rsl(modal) − rsl(VE) [m]",
                 width = 11, height = Relative(0.42), labelsize = 12, ticklabelsize = 11)
        Label(fig[0, 0:4], @sprintf("%s — rsl at %.0f ka  (modal-N: %s)", set, tka, rank);
              fontsize = 19, font = :bold)
        resize_to_layout!(fig)
        tag = @sprintf("%02dka", round(Int, abs(tka)))
        fn = joinpath(OUTDIR, "rsl_panels_$(set)_$(tag).png")
        save(fn, fig); println("wrote ", fn)
    end
end

# ---------------------------------------------------------------------------
function main()
    mkpath(OUTDIR)
    R = load(JLD2_IN, "results")
    sets = R["sets"]
    fig_convergence(R, sets)
    fig_pareto(R, sets)
    fig_bsl(R, sets)
    fig_resid_maps(R, sets)
    fig_rsl_panels(R, RSL_SET)
end

main()
