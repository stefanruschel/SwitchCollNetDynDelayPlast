# export_figure1.jl
#
# Turns the simulation output of simulate_figure1.jl into the plain-text tables that Figure1.tex (pgfplots) draws.
# No simulation is done here, so changes to the appearance of the figure never require rerunning a simulation.
#
# Usage (from this folder):
#   julia --project=.. export_figure1.jl        then        lualatex Figure1.tex
#
# Inputs:  data/figure1b.csv, data/figure1c.csv, data/figure1d.csv   (simulate_figure1.jl)
#          ../connectome_data/sc_raw.mat, D.mat                      (column (a))
# Outputs: tables/*.dat and tables/params.tex (time windows, axis ranges and ticks as TeX macros)
#
# What is drawn:
#   column (a)   (a1) weights, (a2) distances (grey scale, normalised to [0, 1]), (a3) weights against raw distances
#                with an exponential least-squares fit
#   columns 1, 2 raster plots: a tick whenever the phase of a node crosses a multiple of 2 pi,
#                column 1 for t in [0, 100], column 2 for the 100 time units given in SECOND_WINDOW
#   column 3     (b3) mean frequencies over [1500, 2000]; (c3), (d3) mean speed against collective frequency for t >= 1150
#   column 4     all speeds (every 4 time units) and the order parameters R1, R2 (every time unit)
#
# Order of the nodes in the raster plots:
#   (b1), (c1)   by index
#   (b2), (b3)   by mean frequency over [1500, 2000], fastest in the middle (bell_order)
#   (c2), (d1)   by descending phase mod 2 pi at t = 100
#   (d2)         by descending phase difference to the mean phase at t = 1240 (middle of its window)

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
using DelimitedFiles, Printf

const TABLES = joinpath(@__DIR__, "tables")
const FIRST_WINDOW = (0.0, 100.0)
const SECOND_WINDOW = Dict("b" => (1900.0, 2000.0),     # chimera state
                           "c" => (1110.0, 1210.0),     # incoherent phase and resynchronisation
                           "d" => (1190.0, 1290.0))     # incoherent phase and formation of a new two-cluster state

# ------------------------------------------------------------------------------------------------ helpers

"Reads a CSV file written by simulate_figure1.jl."
function load_run(file)
    table, _ = readdlm(file, ','; header = true)
    N = (size(table, 2) - 1) ÷ 2
    return (t = table[:, 1], phases = permutedims(table[:, 2:N+1]), speeds = permutedims(table[:, N+2:end]), N = N)
end

"Times in [t0, t1] at which the unwrapped phase of node i crosses a multiple of 2 pi (linear interpolation)."
function spike_times(run, i, t0, t1)
    t = run.t
    phase = view(run.phases, i, :)
    spikes = Float64[]
    for n in 1:length(t)-1
        if t[n+1] < t0 || t[n] > t1
            continue
        end
        for k in (floor(phase[n] / 2pi) + 1):floor(phase[n+1] / 2pi)
            crossing = t[n] + (2pi * k - phase[n]) / (phase[n+1] - phase[n]) * (t[n+1] - t[n])
            if t0 <= crossing <= t1
                push!(spikes, crossing)
            end
        end
    end
    return spikes
end

"Index of the saved time closest to t."
index_at(run, t) = argmin(abs.(run.t .- t))

"Mean frequency of every node over [t0, t1]."
function mean_frequency(run, t0, t1)
    i0 = index_at(run, t0)
    i1 = index_at(run, t1)
    return (run.phases[:, i1] .- run.phases[:, i0]) ./ (run.t[i1] - run.t[i0])
end

"Nodes ordered by descending phase mod 2 pi at time t."
order_by_phase(run, t) = sortperm(mod.(run.phases[:, index_at(run, t)], 2pi); rev = true)

"Nodes ordered by descending phase difference to the mean phase arg(sum_j exp(i theta_j)) at time t, wrapped to (-pi, pi]."
function order_by_phase_difference(run, t)
    phase = run.phases[:, index_at(run, t)]
    mean_phase = angle(sum(exp.(im .* phase)))
    return sortperm(angle.(exp.(im .* (phase .- mean_phase))); rev = true)
end

"Writes a table with a header line, numbers rounded to 5 significant digits."
function write_table(name, header, matrix)
    open(joinpath(TABLES, name), "w") do io
        println(io, header)
        writedlm(io, round.(matrix; sigdigits = 5), ' ')
    end
end

"""
Mean speed and collective frequency (central difference of the mean phase over +-h) for t >= tmin, every `every` time units:
the projections of Fig. 1(c3), (d3).
"""
function speed_frequency_projection(run; tmin = 1150.0, h = 2.0, every = 1.0)
    mean_phase = vec(mean(run.phases, dims = 1))
    mean_speed = vec(mean(run.speeds, dims = 1))
    dt = run.t[2] - run.t[1]
    shift = round(Int, h / dt)
    stride = max(1, round(Int, every / dt))
    indices = [i for i in shift+1:stride:length(run.t)-shift if run.t[i] >= tmin]
    frequency = [(mean_phase[i+shift] - mean_phase[i-shift]) / (2shift * dt) for i in indices]
    return hcat(mean_speed[indices], frequency)
end

"Comma-separated tick positions with a round step, about three ticks in [lo, hi]."
function nice_ticks(lo, hi)
    steps = [0.05, 0.1, 0.2, 0.25, 0.5, 1.0]
    step = steps[findfirst(s -> (hi - lo) / s <= 3.5, steps)]
    ticks = collect(ceil(lo / step) * step:step:hi)
    return join(string.(round.(ticks; digits = 2)), ",")
end

"Interval [lo, hi] enlarged by a fraction f of its length on both sides."
padded(lo, hi; f = 0.06) = (lo - f * (hi - lo), hi + f * (hi - lo))

# ------------------------------------------------------------------------------------------------ export

function main()
    mkpath(TABLES)
    runs = Dict(row => load_run(joinpath(@__DIR__, "data", "figure1$row.csv")) for row in ("b", "c", "d"))
    N = runs["b"].N

    frequencies_b = mean_frequency(runs["b"], 1500.0, 2000.0)
    node_order = Dict(("b", 1) => collect(1:N), ("b", 2) => bell_order(frequencies_b),
                      ("c", 1) => collect(1:N), ("c", 2) => order_by_phase(runs["c"], 100.0),
                      ("d", 1) => order_by_phase(runs["d"], 100.0), ("d", 2) => order_by_phase_difference(runs["d"], 1240.0))

    # columns 1 and 2: raster plots, rows (t, position) with position N at the top and 1 at the bottom
    for row in ("b", "c", "d"), column in (1, 2)
        run = runs[row]
        order = node_order[(row, column)]
        window = column == 1 ? FIRST_WINDOW : SECOND_WINDOW[row]
        points = Float64[]
        for (position, node) in enumerate(order)
            for t in spike_times(run, node, window...)
                append!(points, (t, N + 1 - position))
            end
        end
        write_table("raster_$(row)$(column).dat", "t y", permutedims(reshape(points, 2, :)))
    end

    # column 3: (b3) mean frequencies in the order of (b2); (c3), (d3) projections onto (mean speed, frequency)
    order = node_order[("b", 2)]
    write_table("freq_b3.dat", "omega y", hcat(frequencies_b[order], N + 1 .- (1:N)))
    projection_c = speed_frequency_projection(runs["c"])
    projection_d = speed_frequency_projection(runs["d"])
    write_table("traj_c.dat", "v omega", projection_c)
    write_table("traj_d.dat", "v omega", projection_d)

    # column 4: speeds every 4 time units, order parameters every time unit
    for row in ("b", "c", "d")
        run = runs[row]
        dt = run.t[2] - run.t[1]
        speed_indices = 1:round(Int, 4 / dt):length(run.t)
        order_indices = 1:round(Int, 1 / dt):length(run.t)
        speed_header = join(vcat("t", ["v$i" for i in 1:N]), ' ')
        write_table("speeds_$(row).dat", speed_header, hcat(run.t[speed_indices], permutedims(run.speeds[:, speed_indices])))
        R1 = order_parameter(run.phases[:, order_indices], 1)
        R2 = order_parameter(run.phases[:, order_indices], 2)
        write_table("order_$(row).dat", "t R1 R2", hcat(run.t[order_indices], R1, R2))
    end

    # column (a): connectome, weights with average row sum 1 and distances with mean 1
    raw_weights, raw_distances = raw_connectome()
    weights = raw_weights ./ mean(sum(raw_weights, dims = 2))
    distances = raw_distances ./ mean(raw_distances)
    for (name, matrix) in (("weights.dat", weights ./ maximum(weights)), ("distances.dat", distances ./ maximum(distances)))
        entries = reduce(vcat, [[j i matrix[i, j]] for i in 1:N for j in 1:N])
        write_table(name, "j i z", entries)
    end
    # (a3) weights against raw distances (off-diagonal pairs with positive weight), fit w = exp(p1 d + p2)
    pairs = (raw_weights .> 0) .& (raw_distances .> 0) .& .!(Matrix(I, N, N))
    x = raw_distances[pairs]
    y = weights[pairs]
    fit = hcat(x, ones(length(x))) \ log.(y)
    x_fit = range(minimum(x), maximum(x); length = 60)
    write_table("scatter_a3.dat", "d w", hcat(x, y))
    write_table("fit_a3.dat", "d w", hcat(x_fit, exp.(fit[1] .* x_fit .+ fit[2])))
    println("exponential fit of the weights against raw distances: exponent $(round(fit[1], digits = 2))")

    # time windows, axis ranges and ticks for Figure1.tex
    open(joinpath(TABLES, "params.tex"), "w") do io
        for row in ("b", "c", "d")
            K = uppercase(row)
            window = SECOND_WINDOW[row]
            @printf(io, "\\def\\WinTwoStart%s{%g}\\def\\WinTwoEnd%s{%g}\n", K, window[1], K, window[2])
            @printf(io, "\\def\\Vmax%s{%.2f}\n", K, 1.06 * maximum(runs[row].speeds))
        end
        for (K, projection) in (("C", projection_c), ("D", projection_d))
            (vlo, vhi) = padded(extrema(projection[:, 1])...)
            (olo, ohi) = padded(extrema(projection[:, 2])...)
            @printf(io, "\\def\\TrajVmin%s{%.3f}\\def\\TrajVmax%s{%.3f}\\def\\TrajOmin%s{%.3f}\\def\\TrajOmax%s{%.3f}\n",
                    K, vlo, K, vhi, K, olo, K, ohi)
            println(io, "\\def\\TrajVticks$(K){$(nice_ticks(vlo, vhi))}\\def\\TrajOticks$(K){$(nice_ticks(olo, ohi))}")
        end
        (flo, fhi) = padded(extrema(frequencies_b)...)
        @printf(io, "\\def\\FreqMin{%.3f}\\def\\FreqMax{%.3f}\n", flo, fhi)
        println(io, "\\def\\FreqTicks{$(nice_ticks(flo, fhi))}")
        @printf(io, "\\def\\DistMax{%.3f}\n", 1.03 * maximum(x))
        @printf(io, "\\def\\WMax{%.4f}\n", 1.03 * maximum(y))
    end
    println("wrote tables/*.dat and tables/params.tex")
end

main()
