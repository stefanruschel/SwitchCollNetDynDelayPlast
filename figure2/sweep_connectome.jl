# sweep_connectome.jl
#
# Functions for the parameter sweeps of Fig. 2: long-time behaviour of the connectome model over a grid of
# coupling strengths sigma and thresholds nu, for a given interaction function. Included by run_figure2_sweeps.jl;
# not meant to be run on its own.
#
# Every run starts from a random initial condition drawn from MersenneTwister(seed) (theta_i ~ U(0, 2 pi),
# v_i ~ N(2.5, 0.1)); seed 1 gives the same initial condition in every grid cell and for every interaction function.
# Adaptation (sigmoid rule (5)) is on from t = 0. All statistics are taken over the second half of the run.
#
# Output: a CSV file with one line per run and the columns
#   sigma, nu, seed       parameters of the run
#   R1mean, R2mean        time averages <R1>, <R2> of the order parameters
#   R1min, R1max, R2min, R2max
#   freq_std              standard deviation of the mean frequencies of the nodes
#   events_R1, events_R2  number of collapses of R1, R2 (rise above 0.8, then fall below 0.5)  [not used in the Letter]
#   vbar_end, vmin_end, vmax_end   mean, smallest and largest speed at the end of the run
#   Omegabar_end          collective frequency at the end of the run (over the last 10 time units)
#   wall_s                computing time of the run in seconds (differs from machine to machine)
#   frac_locked           fraction of nodes whose mean frequency is within 1e-3 of the median mean frequency
#   events_R1_strict, events_R2_strict   number of collapses of R1, R2 (rise above 0.9, then fall below 0.1);
#                         these define "switching" in Fig. 2 (at least two collapses of R1 or R2)

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
using DelimitedFiles, Printf

# Solver tolerance used for the published sweeps of Fig. 2 (abstol = reltol). See the README: at this tolerance
# single nodes can make spurious 2 pi phase slips, which is why the time series of Fig. 1 were computed at 1e-7.
const SWEEP_TOLERANCE = 1e-4

# The six interaction functions of Fig. 2, h(x) = amplitude (-sin(x - a) + r sin 2x), with the parameter values
# exactly as used for the published runs (a = pi + phase lag, rounded to 8 decimals).
const INTERACTION_FUNCTIONS = Dict(
    "b" => (formula = "sin x", a = 3.14159265, r = 0.0, amplitude = 1.0),
    "c" => (formula = "sin(x - pi/4)", a = 3.92699082, r = 0.0, amplitude = 1.0),
    "d" => (formula = "sin(x - 3pi/8)", a = 4.31968990, r = 0.0, amplitude = 1.0),
    "e" => (formula = "sin(x - pi/2)", a = 4.71238898, r = 0.0, amplitude = 1.0),
    "f" => (formula = "c_1 [sin(x - pi/2) + sin(2x)/3], c_1 = 0.86083339", a = 4.71238898, r = 0.33333333, amplitude = 0.86083339),
    "g" => (formula = "c_2 [sin(x - pi/2) + sin 2x], c_2 = 0.56812611", a = 4.71238898, r = 1.0, amplitude = 0.56812611),
)

const CSV_HEADER = "sigma,nu,seed,R1mean,R2mean,R1min,R1max,R2min,R2max,freq_std,events_R1,events_R2," *
                   "vbar_end,vmin_end,vmax_end,Omegabar_end,wall_s,frac_locked,events_R1_strict,events_R2_strict"

"""
    simulate_and_summarise(W, D, sigma, nu, seed; a, r, amplitude, tend)

One run of the connectome model; returns the numbers of one line of the sweep CSV (see the header of this file).
"""
function simulate_and_summarise(W, D, sigma, nu, seed; a, r, amplitude, tend)
    p = ConnectomeParameters(W = W, D = D, sigma = sigma, nu = nu, a = a, r = r, amplitude = amplitude)
    N = size(W, 1)
    rng = MersenneTwister(seed)
    theta0, v0 = random_initial_condition(rng, N; vmin = p.vmin, vmax = p.vmax)

    wall_time = @elapsed solution = simulate(p, theta0, v0; tend = tend, saveat = 0.5, tolerance = SWEEP_TOLERANCE)

    times, phases, speeds = solution_window(solution, N, tend / 2, tend)      # second half of the run
    R1 = order_parameter(phases, 1)
    R2 = order_parameter(phases, 2)
    frequencies = mean_frequencies(times, phases)
    collective_frequency_end = (mean(phases[:, end]) - mean(phases[:, end-20])) / (times[end] - times[end-20])

    return [sigma, nu, seed,
            mean(R1), mean(R2), minimum(R1), maximum(R1), minimum(R2), maximum(R2),
            std(frequencies),
            count_collapses(R1; high = 0.8, low = 0.5), count_collapses(R2; high = 0.8, low = 0.5),
            mean(speeds[:, end]), minimum(speeds[:, end]), maximum(speeds[:, end]),
            collective_frequency_end, wall_time, fraction_locked(frequencies),
            count_collapses(R1; high = 0.9, low = 0.1), count_collapses(R2; high = 0.9, low = 0.1)]
end

"""
    run_sweep(output_file; interaction, sigma_min, sigma_max, n_sigma, nu_min, nu_max, n_nu, n_seeds, tend)

Runs the model on an n_sigma x n_nu grid (sigma logarithmically spaced in [sigma_min, sigma_max], nu linearly
spaced in [nu_min, nu_max]) with seeds 1..n_seeds per grid cell, using all Julia threads, and writes the CSV file.
`interaction` is one entry of INTERACTION_FUNCTIONS.
"""
function run_sweep(output_file; interaction, sigma_min, sigma_max, n_sigma, nu_min, nu_max, n_nu, n_seeds, tend)
    sigmas = exp.(range(log(sigma_min), log(sigma_max); length = n_sigma))
    nus = range(nu_min, nu_max; length = n_nu)
    W, D = load_connectome()

    jobs = vec([(sigma, nu, seed) for sigma in sigmas, nu in nus, seed in 1:n_seeds])
    println("h(x) = $(interaction.formula): $(length(jobs)) runs on $(Threads.nthreads()) threads, t <= $tend -> $output_file")
    results = Vector{Vector{Float64}}(undef, length(jobs))
    Threads.@threads for k in eachindex(jobs)
        sigma, nu, seed = jobs[k]
        results[k] = simulate_and_summarise(W, D, sigma, nu, seed;
                                            a = interaction.a, r = interaction.r, amplitude = interaction.amplitude, tend = tend)
        @printf("  sigma = %.4f  nu = %.3f  seed = %d   <R1> = %.2f  <R2> = %.2f  locked = %.2f  collapses (R1, R2) = (%d, %d)  (%.0f s)\n",
                sigma, nu, seed, results[k][4], results[k][5], results[k][18], results[k][19], results[k][20], results[k][17])
    end

    open(output_file, "w") do io
        println(io, CSV_HEADER)
        writedlm(io, reduce(vcat, (line' for line in results)), ',')
    end
    println("wrote $output_file")
end
