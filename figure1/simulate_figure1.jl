# simulate_figure1.jl
#
# Simulations behind Fig. 1(b)-(d) of the Letter (68-node connectome, sigmoid adaptation rule (5)).
#
# Usage (from this folder):
#   julia --project=.. simulate_figure1.jl b        row (b): chimera state
#   julia --project=.. simulate_figure1.jl c        row (c): switching between synchrony and incoherence
#   julia --project=.. simulate_figure1.jl d        row (d): switching between two-cluster states
#   julia --project=.. simulate_figure1.jl b --quick     short test run (t <= 60), written to data/quick_figure1b.csv
#
# Run time: about one minute per row on one core (plus a few minutes of compilation on the first run).
#
# Output: data/figure1<row>.csv with the header
#   t, theta_1, ..., theta_68, v_1, ..., v_68
# one line every 0.5 time units for 0 <= t <= 2000 (phases unwrapped, i.e. not taken mod 2 pi),
# and data/figure1<row>.params.txt with the parameters of the run.

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))
using DelimitedFiles

# Parameters that differ between the rows of Fig. 1. The values of a and amplitude are written exactly as used for
# the published runs (a = pi + phase lag; see src/connectome_model.jl):
#   (b), (c):  h(x) = sin x                                   a = pi
#   (d):       h(x) = c_2 [sin(x - pi/2) + sin 2x]            a = 4.71238898 (= 3 pi / 2), r = 1, c_2 = 0.56812611
const ROWS = Dict(
    "b" => (description = "chimera state",
            sigma = 0.6, nu = 1.44, a = pi, r = 0.0, amplitude = 1.0,
            t_on = 1000.0),                       # speeds kept constant until t = 1000
    "c" => (description = "switching between synchrony and incoherence",
            sigma = 0.8, nu = 1.56, a = pi, r = 0.0, amplitude = 1.0,
            t_on = 0.0),
    "d" => (description = "switching between two-cluster states",
            sigma = 1.0, nu = 1.4, a = 4.71238898, r = 1.0, amplitude = 0.56812611,
            t_on = 0.0),
)

# Parameters common to all rows
const SEED = 1                 # random initial condition: Random.seed!(1), then theta_i ~ U(0, 2 pi), v_i ~ N(2.5, 0.1)
const TOLERANCE = 1e-7         # abstol = reltol of the solver
const SAVE_STEP = 0.5          # output interval

function main()
    if isempty(ARGS) || !haskey(ROWS, ARGS[1])
        error("give the row of Fig. 1 as first argument: b, c or d")
    end
    row = ARGS[1]
    quick = flag("quick")
    tend = quick ? 60.0 : 2000.0
    run = ROWS[row]
    println("Fig. 1($row): $(run.description)")

    W, D = load_connectome()
    p = ConnectomeParameters(W = W, D = D, sigma = run.sigma, nu = run.nu, a = run.a, r = run.r,
                             amplitude = run.amplitude, t_on = run.t_on)

    Random.seed!(SEED)
    theta0, v0 = random_initial_condition(Random.default_rng(), NUMBER_OF_NODES; vmin = p.vmin, vmax = p.vmax)

    elapsed = @elapsed solution = simulate(p, theta0, v0; tend = tend, saveat = SAVE_STEP, tolerance = TOLERANCE)
    println("integrated to t = $tend in $(round(elapsed / 60, digits = 1)) min")
    times, phases, speeds = solution_window(solution, NUMBER_OF_NODES, 0.0, tend)

    mkpath(joinpath(@__DIR__, "data"))
    name = (quick ? "quick_" : "") * "figure1$row"
    csv_file = joinpath(@__DIR__, "data", name * ".csv")
    header = vcat("t", ["theta_$i" for i in 1:NUMBER_OF_NODES], ["v_$i" for i in 1:NUMBER_OF_NODES])
    open(csv_file, "w") do io
        println(io, join(header, ","))
        writedlm(io, hcat(times, phases', speeds'), ',')
    end
    open(joinpath(@__DIR__, "data", name * ".params.txt"), "w") do io
        println(io, "Fig. 1($row), $(run.description)")
        println(io, "sigma=$(p.sigma) nu=$(p.nu) a=$(p.a) r=$(p.r) amplitude=$(p.amplitude) rule=$(p.rule)")
        println(io, "omega=$(p.omega) epsilon=$(p.epsilon) t_on=$(p.t_on) beta=$(p.beta) window=$(p.window) vmin=$(p.vmin) vmax=$(p.vmax)")
        println(io, "seed=$SEED tolerance=$TOLERANCE tend=$tend saveat=$SAVE_STEP")
    end
    println("wrote $csv_file ($(length(times)) rows)")
end

main()
