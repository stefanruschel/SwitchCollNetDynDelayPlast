# chimera_long_run.jl
#
# Persistence of the chimera state of Fig. 1(b) (SM paragraph "Persistence of the chimera state"): the run of
# simulate_figure1.jl b (same parameters and initial condition, adaptation switched on at t = 1000) continued to
# t = 10^6. In each of the 100 consecutive windows of 10^4 time units the following diagnostics are computed over the
# last 200 time units of the window:
#   R1mean, R1min   mean and minimum of the order parameter R1
#   n_locked        number of nodes whose mean frequency lies within 1e-3 of the median mean frequency
#   freq_std        standard deviation of the mean frequencies of the nodes
#   vbar            mean speed at the end of the window
#
# Usage (from this folder):
#   julia --project=.. chimera_long_run.jl                         t <= 1e6, about 80 minutes on one core
#   julia --project=.. chimera_long_run.jl --tend=2e4 --quick      shorter test run (data/quick_chimera_long_run.csv)
#
# Output: data/chimera_long_run.csv with the columns t (end of the window), R1mean, R1min, n_locked, freq_std, vbar.

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))
using DelimitedFiles, Printf

const TOLERANCE = 1e-4          # solver tolerance of the published run
const WINDOW = 1.0e4            # length of a window
const PROBE = 200.0             # diagnostics over the last PROBE time units of every window

function main()
    tend = option("tend", 1.0e6)
    quick = flag("quick")

    W, D = load_connectome()
    p = ConnectomeParameters(W = W, D = D, sigma = 0.6, nu = 1.44, a = pi, r = 0.0, t_on = 1000.0)
    Random.seed!(1)
    theta0, v0 = random_initial_condition(Random.default_rng(), NUMBER_OF_NODES; vmin = p.vmin, vmax = p.vmax)

    # save every 50 time units, and every 0.5 time units during the probe at the end of each window
    save_times = Float64[]
    for window_start in 0.0:WINDOW:(tend - WINDOW)
        append!(save_times, window_start:50.0:(window_start + WINDOW - PROBE))
        append!(save_times, (window_start + WINDOW - PROBE):0.5:(window_start + WINDOW))
    end
    save_times = unique(sort(save_times))

    elapsed = @elapsed solution = simulate(p, theta0, v0; tend = tend, saveat = save_times, tolerance = TOLERANCE)
    println("integrated to t = $tend in $(round(elapsed / 60, digits = 1)) min")

    rows = Vector{Float64}[]
    for window_start in 0.0:WINDOW:(tend - WINDOW)
        window_end = window_start + WINDOW
        times, phases, speeds = solution_window(solution, NUMBER_OF_NODES, window_end - PROBE, window_end)
        R1 = order_parameter(phases, 1)
        frequencies = mean_frequencies(times, phases)
        n_locked = count(abs.(frequencies .- median(frequencies)) .< 1e-3)
        push!(rows, [window_end, mean(R1), minimum(R1), n_locked, std(frequencies), mean(speeds[:, end])])
        @printf("t = %9.0f   <R1> = %.3f  min R1 = %.3f  locked = %2d of 68  freq_std = %.3e  mean speed = %.3f\n", rows[end]...)
    end

    mkpath(joinpath(@__DIR__, "data"))
    file = joinpath(@__DIR__, "data", (quick ? "quick_" : "") * "chimera_long_run.csv")
    open(file, "w") do io
        println(io, "t,R1mean,R1min,n_locked,freq_std,vbar")
        writedlm(io, reduce(vcat, (row' for row in rows)), ',')
    end
    println("wrote $file")
end

main()
