#!/usr/bin/env julia
# Cross-section and geocentre plots for the 3-D viscosity benchmark (setup note by Klemann et al.)
# (tests/bench_visc3d.f90). For one forcing: a 4x3 panel figure, rows u_r, u_θ,
# u_φ, δφ against distance along the section, columns tests A, B, C, one line per
# output epoch; plus the geocentre components against time. The protocol files
# are searched for anywhere under `dir`, so a runme ensemble directory (one
# member per test/forcing) works as is. Missing tests are left blank.
#
# Usage:  julia --project=analysis analysis/plot_bench_visc3d.jl [dir] [heav|ramp]
#         (defaults: runs/bench_visc3d heav)

using CairoMakie
using Printf

const TESTS  = ["A", "B", "C"]
const FIELDS = [(5, L"u_r\ \mathrm{[m]}"), (6, L"u_\vartheta\ \mathrm{[m]}"),
                (7, L"u_\varphi\ \mathrm{[m]}"), (8, L"\delta\phi\ \mathrm{[m^2\,s^{-2}]}")]
const OUTDIR = "analysis/figs"

# Path of the protocol file `name` under `dir`, or nothing.
function find_file(dir, name)
    for (root, _, files) in walkdir(dir)
        name in files && return joinpath(root, name)
    end
    return nothing
end

# Numeric rows of a protocol file (tab-delimited, '#' header lines).
function read_table(path)
    rows = [parse.(Float64, split(strip(l))) for l in eachline(path)
            if !startswith(l, "#") && !isempty(strip(l))]
    return reduce(vcat, permutedims.(rows))
end

function plot_sections(dir, forcing)
    fig = Figure(size = (1200, 1100))
    for (j, t) in enumerate(TESTS)
        path = find_file(dir, "disp_VILMA2_$(t)-i_$(forcing).txt")
        d = isnothing(path) ? nothing : read_table(path)
        for (i, (col, lab)) in enumerate(FIELDS)
            ax = Axis(fig[i, j]; ylabel = j == 1 ? lab : "",
                      xlabel = i == length(FIELDS) ? "distance from structure centre [deg]" : "",
                      title = i == 1 ? "$(t)-i ($(forcing))" : "")
            xlims!(ax, 0, 180)
            isnothing(d) && continue
            for tk in unique(d[:, 4])
                sel = d[:, 4] .== tk
                lines!(ax, d[sel, 3], d[sel, col]; color = log10(tk),
                       colorrange = (-1, 2), colormap = :viridis)
            end
        end
    end
    Colorbar(fig[:, length(TESTS)+1]; colormap = :viridis, limits = (-1, 2),
             label = "time [kyr]", ticks = (-1:2, ["0.1", "1", "10", "100"]))
    return fig
end

function plot_gcm(dir, forcing)
    fig = Figure(size = (1000, 350))
    for (i, c) in enumerate(["u_x", "u_y", "u_z"])
        ax = Axis(fig[1, i]; xlabel = "time [kyr]", ylabel = "$(c) [m]", xscale = log10,
                  title = "geocentre u_CF − u_CM ($(forcing))")
        for t in TESTS
            path = find_file(dir, "gcm_VILMA2_$(t)-i_$(forcing).txt")
            isnothing(path) && continue
            g = read_table(path)
            scatterlines!(ax, g[:, 1], g[:, i+1]; label = "$(t)-i")
        end
        i == 3 && axislegend(ax; position = :lb)
    end
    return fig
end

function main(args)
    dir     = length(args) >= 1 ? args[1] : "runs/bench_visc3d"
    forcing = length(args) >= 2 ? args[2] : "heav"
    mkpath(OUTDIR)
    tag = replace(dir, r"[/\\]" => "_")
    f1 = joinpath(OUTDIR, "bench_visc3d_sections_$(tag)_$(forcing).png")
    f2 = joinpath(OUTDIR, "bench_visc3d_gcm_$(tag)_$(forcing).png")
    save(f1, plot_sections(dir, forcing))
    save(f2, plot_gcm(dir, forcing))
    @printf("wrote %s\n      %s\n", f1, f2)
end

main(ARGS)
