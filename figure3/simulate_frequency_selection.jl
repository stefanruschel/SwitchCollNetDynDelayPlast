# simulate_frequency_selection.jl
#
# Green trajectory of Fig. 3(a): frequency selection on the ring of N = 11 nodes.
# Initial data (caption of Fig. 3): theta_j(t) = Omega t, v_j(t) = c for t <= 0 with Omega = 0.9996 and c = 1/5.1,
# i.e. on the synchronous branch psi_0 just below the threshold nu = 1. Other parameters: defaults of RingParameters
# (sigma = 0.1, epsilon = 1e-5, nu = 1, v_min = 1/6.8, v_max = 5, h(x) = sin(x - pi/2)).
# The trajectory drifts slowly along psi_0 until the speeds reach v_min.
#
# Usage (from this folder):
#   julia --project=.. simulate_frequency_selection.jl            t <= 3e6, about one minute
#   julia --project=.. simulate_frequency_selection.jl --quick    t <= 2000 (test, files data/quick_*)
#
# Outputs: data/frequency_selection_omega_tau_node1.csv   columns 1/v_1, theta_1'          (drawn in Fig. 3(a))
#          data/frequency_selection_omega_tau_mean.csv    columns mean_j 1/v_j, mean_j theta_j'
# (1000 saved times, frequencies by forward differences between them)

include(joinpath(@__DIR__, "..", "src", "ring_model.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))

function main()
    quick = flag("quick")
    s = RingSystem(RingParameters())
    Omega0 = 0.9996
    c0 = 1 / 5.1
    tend = quick ? 2000.0 : 3e6
    save_times = range(0.0, tend; length = quick ? 2000 : 1000)

    elapsed = @elapsed solution = simulate_from_locked_state(s, Omega0, c0; tend = tend, saveat = save_times)
    println("integrated to t = $tend in $(round(elapsed, digits = 1)) s")

    times, phases, speeds = solution_window(solution, s.p.N, 0.0, tend)
    mkpath(joinpath(@__DIR__, "data"))
    write_omega_tau(joinpath(@__DIR__, "data", (quick ? "quick_" : "") * "frequency_selection"), times, phases, speeds)
end

main()
