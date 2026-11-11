# ring_epsilon_sweep.jl
#
# SM Table S1: dependence of the relaxation cycle of the ring (blue trajectory of Fig. 3) on the adaptation rate
# epsilon. Same parameters and initial data as the blue trajectory (Omega = 1.0004, v = 1/5 on the synchronous state
# psi_0), epsilon varied; each run lasts 40/epsilon (several cycles) with 20000 saved times.
#
# At every saved time the ring counts as close to the m-twisted state psi_m if |Z| > 0.95,
# Z = (1/N) sum_j exp(i (theta_{j+1} - theta_j)) (cyclic), with m from arg Z = 2 pi m / N folded to 0..5.
#
# Usage (from this folder):
#   julia -t auto --project=.. ring_epsilon_sweep.jl       about one minute with 9 threads
#
# Output: data/ring_epsilon_sweep.csv with the columns
#   eps                 adaptation rate
#   tend                run length 40/eps
#   frac_on_branch      fraction of the saved times close to some psi_m, after the first departure from psi_0
#                       ("time near a branch" in Table S1)
#   n_departures_psi0   number of departures from psi_0 in the visit sequence
#   sequence            visit sequence: twist numbers m of the states held for at least 0.2/eps, repeated
#                       consecutive entries merged; "transitions" in Table S1 = number of entries - 1
#   wall_s              computing time in seconds (machine dependent)

include(joinpath(@__DIR__, "..", "src", "ring_model.jl"))
using Printf

const EPSILONS = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2]
const CYCLES = 40.0                 # run length = CYCLES / epsilon

"""
    analyse(solution, N, epsilon)

Fraction of time near a twisted state after the first departure from psi_0, number of departures from psi_0, and the
visit sequence (see the header of this file).
"""
function analyse(solution, N, epsilon)
    times = solution.t
    phases = reduce(hcat, [u[1:N] for u in solution.u])
    twist = twist_number(phases; threshold = 0.95)
    n = length(times)

    first_departure = findfirst(k -> twist[k] != 0 && times[k] > 10.0, 1:n)
    after_departure = first_departure === nothing ? (1:n) : (first_departure:n)
    fraction_on_branch = count(>=(0), twist[after_departure]) / length(after_departure)

    sequence = Int[]
    k = 1
    while k <= n
        if twist[k] >= 0
            j = k
            while j < n && twist[j+1] == twist[k]
                j += 1
            end
            held_long_enough = times[j] - times[k] >= 0.2 / epsilon
            if held_long_enough && (isempty(sequence) || sequence[end] != twist[k])
                push!(sequence, twist[k])
            end
            k = j + 1
        else
            k += 1
        end
    end
    departures = count(i -> sequence[i] == 0 && i < length(sequence), eachindex(sequence))
    return fraction_on_branch, departures, sequence
end

function main()
    results = Vector{Any}(undef, length(EPSILONS))
    Threads.@threads for i in eachindex(EPSILONS)
        epsilon = EPSILONS[i]
        s = RingSystem(RingParameters(epsilon = epsilon))
        tend = CYCLES / epsilon
        wall_time = @elapsed solution = simulate_from_locked_state(s, 1.0004, 1 / 5; tend = tend,
                                                                   saveat = range(0.0, tend; length = 20000))
        fraction, departures, sequence = analyse(solution, s.p.N, epsilon)
        results[i] = (epsilon, tend, fraction, departures, sequence, wall_time)
        @printf("eps = %.0e  t <= %.3g  time near a branch %.3f  departures from psi_0 %d  transitions %d  sequence %s  (%.0f s)\n",
                epsilon, tend, fraction, departures, length(sequence) - 1, join(sequence, "-"), wall_time)
    end
    mkpath(joinpath(@__DIR__, "data"))
    file = joinpath(@__DIR__, "data", "ring_epsilon_sweep.csv")
    open(file, "w") do io
        println(io, "eps,tend,frac_on_branch,n_departures_psi0,sequence,wall_s")
        for r in results
            @printf(io, "%.3g,%.4g,%.4f,%d,%s,%.0f\n", r[1], r[2], r[3], r[4], join(r[5], "-"), r[6])
        end
    end
    println("wrote $file")
end

main()
