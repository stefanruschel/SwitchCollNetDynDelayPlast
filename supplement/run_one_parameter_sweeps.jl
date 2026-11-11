# run_one_parameter_sweeps.jl
#
# One-parameter sweeps about the parameters of Fig. 1(d) for SM Fig. S1(b)-(f) ("One-parameter sweeps" in the SM).
# Base point: sigma = 1.0, nu = 1.4, h(x) = c_2 [sin(x - pi/2) + sin 2x], sigmoid rule (5) with beta = 5,
# epsilon = 0.01, v_i(0) ~ N(2.5, 0.1), weights and distances as in Fig. 1. Six initial conditions (seeds 1..6,
# MersenneTwister) per parameter value. Each run lasts until t = max(2000, K/epsilon) (K = 0 or 20, see SWEEPS),
# and all statistics are taken over the second half of the run.
#
# Usage (from this folder; -t auto uses all cores):
#   julia -t auto --project=.. run_one_parameter_sweeps.jl                 all nine sweeps, about one hour on 13 threads
#   julia -t auto --project=.. run_one_parameter_sweeps.jl --only=eps,beta  some of them (names as in SWEEPS)
#
# Output: data/sweep1d_<name>.csv with one line per run and the columns
#   value               value of the swept parameter
#   seed                initial condition
#   R1mean, R2mean      time averages <R1>, <R2>
#   R2min               minimum of R2
#   freq_std            standard deviation of the mean frequencies of the nodes
#   events_R2           collapses of R2 (rise above 0.8, then fall below 0.5)  [not used in the SM]
#   vbar_end            mean speed at the end of the run
#   wall_s              computing time in seconds (machine dependent)
#   tend                length of the run
#   events_R1s, events_R2s   collapses of R1, R2 (rise above 0.9, then fall below 0.1); a run "switches" if one of
#                            them is at least 2 (as in Fig. 2)

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))
using DelimitedFiles, Printf

const SOLVER_TOLERANCE = 1e-4          # as for the sweeps of Fig. 2
const NUMBER_OF_SEEDS = 6

# base point = Fig. 1(d)
const BASE = (sigma = 1.0, nu = 1.4, a = 4.71238898, r = 1.0, amplitude = 0.56812611,
              beta = 5.0, epsilon = 0.01, speed_spread = 0.1, rule = :sigmoid)

# The sweeps. `parameter` is the swept quantity:
#   :epsilon        adaptation rate                             (Fig. S1(b))
#   :speed_spread   standard deviation s of v_i(0) ~ N(2.5, s)  (Fig. S1(c))
#   :beta           gain of the sigmoid                         (Fig. S1(d))
#   :nu             threshold of rule (5), other beta           (Fig. S1(e))
#   :rule6_nu       threshold of the polynomial rule (6)        (Fig. S1(f))
# `changes` modifies the base point; `K` sets the run length max(2000, K/epsilon).
const SWEEPS = [
    (name = "eps", parameter = :epsilon, values = [0.001, 0.003, 0.01, 0.03], changes = (;), K = 20.0),
    (name = "spread", parameter = :speed_spread, values = [0.01, 0.1, 0.5, 1, 2], changes = (;), K = 0.0),
    (name = "beta", parameter = :beta, values = [1, 2, 3, 5, 8, 12], changes = (;), K = 0.0),
    (name = "nu_beta3", parameter = :nu, values = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2],
     changes = (beta = 3.0,), K = 0.0),
    (name = "nu_beta8", parameter = :nu, values = [0.9, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5],
     changes = (beta = 8.0,), K = 0.0),
    (name = "rule6_nu", parameter = :rule6_nu, values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6],
     changes = (;), K = 0.0),
    (name = "rule6_nu_fine", parameter = :rule6_nu, values = [0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98],
     changes = (;), K = 0.0),
    (name = "rule6_nu_eps0.002", parameter = :rule6_nu, values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
     changes = (epsilon = 0.002,), K = 20.0),
    (name = "rule6_nu_eps0.002_fine", parameter = :rule6_nu, values = [0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98],
     changes = (epsilon = 0.002,), K = 20.0),
]

const CSV_HEADER = "value,seed,R1mean,R2mean,R2min,freq_std,events_R2,vbar_end,wall_s,tend,events_R1s,events_R2s"

"Parameters of one run: the base point with the sweep's changes and the swept parameter set to `value`."
function run_settings(sweep, value)
    settings = merge(BASE, sweep.changes)
    if sweep.parameter == :rule6_nu
        settings = merge(settings, (rule = :polynomial, nu = Float64(value)))
    else
        settings = merge(settings, NamedTuple{(sweep.parameter,)}((Float64(value),)))
    end
    return settings
end

"One run; returns the numbers of one line of the CSV file."
function simulate_and_summarise(W, D, sweep, value, seed)
    settings = run_settings(sweep, value)
    tend = max(2000.0, sweep.K / settings.epsilon)
    p = ConnectomeParameters(W = W, D = D, sigma = settings.sigma, nu = settings.nu, a = settings.a, r = settings.r,
                             amplitude = settings.amplitude, rule = settings.rule, epsilon = settings.epsilon,
                             beta = settings.beta)
    rng = MersenneTwister(seed)
    theta0, v0 = random_initial_condition(rng, size(W, 1); vmin = p.vmin, vmax = p.vmax, speed_spread = settings.speed_spread)

    wall_time = @elapsed solution = simulate(p, theta0, v0; tend = tend, saveat = 0.5, tolerance = SOLVER_TOLERANCE)

    times, phases, speeds = solution_window(solution, size(W, 1), tend / 2, tend)
    R1 = order_parameter(phases, 1)
    R2 = order_parameter(phases, 2)
    frequencies = mean_frequencies(times, phases)
    return [value, seed, mean(R1), mean(R2), minimum(R2), std(frequencies),
            count_collapses(R2; high = 0.8, low = 0.5), mean(speeds[:, end]), wall_time, tend,
            count_collapses(R1; high = 0.9, low = 0.1), count_collapses(R2; high = 0.9, low = 0.1)]
end

function run_one_parameter_sweep(sweep, W, D)
    output_file = joinpath(@__DIR__, "data", "sweep1d_$(sweep.name).csv")
    jobs = vec([(value, seed) for value in sweep.values, seed in 1:NUMBER_OF_SEEDS])
    println("sweep $(sweep.name): $(length(jobs)) runs on $(Threads.nthreads()) threads -> $output_file")
    results = Vector{Vector{Float64}}(undef, length(jobs))
    Threads.@threads for k in eachindex(jobs)
        value, seed = jobs[k]
        results[k] = simulate_and_summarise(W, D, sweep, value, seed)
        @printf("  %s = %.4g  seed = %d   <R1> = %.2f  <R2> = %.2f  collapses (R1, R2) = (%d, %d)  t <= %.0f  (%.0f s)\n",
                sweep.name, value, seed, results[k][3], results[k][4], results[k][11], results[k][12], results[k][10], results[k][9])
    end
    open(output_file, "w") do io
        println(io, CSV_HEADER)
        writedlm(io, reduce(vcat, (line' for line in results)), ',')
    end
    println("wrote $output_file")
end

function main()
    mkpath(joinpath(@__DIR__, "data"))
    W, D = load_connectome()
    selected = split(option("only", join([s.name for s in SWEEPS], ",")), ",")
    for sweep in SWEEPS
        if sweep.name in selected
            run_one_parameter_sweep(sweep, W, D)
        end
    end
end

main()
