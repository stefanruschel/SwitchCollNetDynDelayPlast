# simulate_figureS1a.jl
#
# Time series of SM Fig. S1(a): the polynomial adaptation rule (6) instead of the sigmoid rule (5) at the parameters of
# Fig. 1(d), with threshold nu = 0.9 and slower adaptation epsilon = 0.002; h(x) = c_2 [sin(x - pi/2) + sin 2x].
# Same initial condition as Fig. 1 (Random.seed!(1), theta_i ~ U(0, 2 pi), v_i ~ N(2.5, 0.1)).
#
# Usage (from this folder):
#   julia --project=.. simulate_figureS1a.jl            t <= 10000, a few minutes
#   julia --project=.. simulate_figureS1a.jl --quick    t <= 60 (test, data/quick_figureS1a.csv)
#
# Output: data/figureS1a.csv with the header t, theta_1, ..., theta_68, v_1, ..., v_68 (one line per time unit).

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))
using DelimitedFiles

const TOLERANCE = 1e-4      # solver tolerance of the published run (as for the sweeps)

function main()
    quick = flag("quick")
    tend = quick ? 60.0 : 10000.0
    W, D = load_connectome()
    p = ConnectomeParameters(W = W, D = D, sigma = 1.0, nu = 0.9, a = 4.71238898, r = 1.0, amplitude = 0.56812611,
                             rule = :polynomial, epsilon = 0.002)
    Random.seed!(1)
    theta0, v0 = random_initial_condition(Random.default_rng(), NUMBER_OF_NODES; vmin = p.vmin, vmax = p.vmax)

    elapsed = @elapsed solution = simulate(p, theta0, v0; tend = tend, saveat = 1.0, tolerance = TOLERANCE)
    println("integrated to t = $tend in $(round(elapsed / 60, digits = 1)) min")
    times, phases, speeds = solution_window(solution, NUMBER_OF_NODES, 0.0, tend)

    mkpath(joinpath(@__DIR__, "data"))
    file = joinpath(@__DIR__, "data", (quick ? "quick_" : "") * "figureS1a.csv")
    header = vcat("t", ["theta_$i" for i in 1:NUMBER_OF_NODES], ["v_$i" for i in 1:NUMBER_OF_NODES])
    open(file, "w") do io
        println(io, join(header, ","))
        writedlm(io, hcat(times, phases', speeds'), ',')
    end
    println("wrote $file ($(length(times)) rows)")
end

main()
