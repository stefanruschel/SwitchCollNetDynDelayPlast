# ring_model.jl
#
# The model of the Letter on a non-local ring of N = 11 nodes (Fig. 3 and SM Table S1):
#
#   theta_i'(t) = omega + sigma * sum_j w_ij h( theta_j(t - tau_ij) - theta_i(t) ),   tau_ij = d_ij / v_j,
#   v_j'(t)     = epsilon (theta_j' - nu) (v_j - v_min) (v_max - v_j)                  (polynomial rule, Eq. (6)),
#
# with ring distances d_ij = min(|i - j|, N - |i - j|), weights w_ij = exp(-d_ij) without self-coupling and
# normalised to unit row sum, and the interaction function h(x) = -sin(x - a) + r sin(2x).
# The default a = 3 pi / 2, r = 0 gives h(x) = -sin(x - 3 pi / 2) = sin(x - pi / 2), as in the caption of Fig. 3.
#
# Numerical method: DelayDiffEq.jl, MethodOfSteps(Tsit5()), abstol = reltol = 1e-8.
# The history for t <= 0 is a phase-locked state theta_j(t) = phi_j + Omega t with constant speeds v_j = c.

using DifferentialEquations, DelayDiffEq, LinearAlgebra, Statistics, Random, DelimitedFiles

"Interaction function h(x) = -sin(x - a) + r sin(2x)."
interaction_function(x, a, r) = -sin(x - a) + r * sin(2x)

"Parameters of the ring (defaults = caption of Fig. 3)."
Base.@kwdef struct RingParameters
    N::Int = 11              # number of nodes
    omega::Float64 = 1.0     # natural frequency
    sigma::Float64 = 0.1     # coupling strength
    a::Float64 = 3pi / 2     # phase parameter of h
    r::Float64 = 0.0         # weight of the second harmonic of h
    epsilon::Float64 = 1e-5  # adaptation rate
    nu::Float64 = 1.0        # threshold
    vmin::Float64 = 1 / 6.8  # lower bound of the speeds
    vmax::Float64 = 5.0      # upper bound of the speeds
end

"Ring distances d_ij = min(|i - j|, N - |i - j|)."
ring_distances(N) = [min(abs(i - j), N - abs(i - j)) for i in 1:N, j in 1:N]

"Weights w_ij = exp(-d_ij), w_ii = 0, normalised to unit row sum."
function ring_weights(distances)
    weights = exp.(-Float64.(distances))
    weights[diagind(weights)] .= 0.0
    return weights ./ sum(weights[1, :])
end

"Parameters together with the coupling matrices and the largest possible delay."
struct RingSystem
    p::RingParameters
    W::Matrix{Float64}
    D::Matrix{Float64}
    maximum_delay::Float64
end

function RingSystem(p::RingParameters)
    distances = ring_distances(p.N)
    return RingSystem(p, ring_weights(distances), Float64.(distances), maximum(distances) / p.vmin)
end

"Right-hand side of the delay differential equation; state u = (theta_1..theta_N, v_1..v_N)."
function ring_rhs!(du, u, history, s::RingSystem, t)
    p = s.p
    N = p.N
    theta = view(u, 1:N)
    speed = view(u, N+1:2N)
    dtheta = view(du, 1:N)
    dspeed = view(du, N+1:2N)
    for i in 1:N
        coupling = 0.0
        for j in 1:N
            if s.W[i, j] != 0.0
                delay = s.D[i, j] / max(speed[j], p.vmin)
                delayed_phase = history(s, t - delay)[j]
                coupling += s.W[i, j] * interaction_function(delayed_phase - theta[i], p.a, p.r)
            end
        end
        dtheta[i] = p.omega + p.sigma * coupling
    end
    for j in 1:N
        dspeed[j] = p.epsilon * (dtheta[j] - p.nu) * (speed[j] - p.vmin) * (p.vmax - speed[j])
    end
    return nothing
end

"""
    simulate_from_locked_state(s, Omega, c; tend, saveat, phase_offsets = zeros(N))

Integrates the ring from the history theta_j(t) = phase_offsets_j + Omega t, v_j(t) = c for t <= 0
(phase_offsets = 0 is the synchronous state psi_0).
"""
function simulate_from_locked_state(s::RingSystem, Omega, c; tend, saveat, phase_offsets = zeros(s.p.N))
    N = s.p.N
    history = (system, t) -> vcat(phase_offsets .+ Omega .* t, fill(c, N))
    problem = DDEProblem(ring_rhs!, history(s, 0.0), history, (0.0, tend), s; constant_lags = [1.0, s.maximum_delay])
    return solve(problem, MethodOfSteps(Tsit5()); abstol = 1e-8, reltol = 1e-8, saveat = saveat)
end

# ------------------------------------------------------------------------------------------------ diagnostics

"Saved times in [tmin, tmax] and the corresponding phases and speeds (N x n matrices)."
function solution_window(solution, N, tmin, tmax)
    indices = findall(t -> tmin <= t <= tmax, solution.t)
    states = reduce(hcat, solution.u[indices])
    return solution.t[indices], states[1:N, :], states[N+1:2N, :]
end

"Kuramoto order parameter R(t) = |sum_j exp(i theta_j)| / N for every column (time) of `phases`."
order_parameter(phases) = [abs(sum(exp.(im .* phases[:, k]))) / size(phases, 1) for k in 1:size(phases, 2)]

"""
Mean nearest-neighbour coherence C(t) = (1/N) sum_i (1 + cos(theta_{i+1} - theta_i)) / 2, cyclic (theta_{N+1} = theta_1),
for every column (time) of `phases`; 0 <= C <= 1 (Fig. 3(c)).
"""
nearest_neighbour_coherence(phases) =
    [sum(1.0 .+ cos.(phases[:, k] .- circshift(phases[:, k], -1))) / (2size(phases, 1)) for k in 1:size(phases, 2)]

"""
    twist_number(phases; threshold = 0.95)

For every column (time) of `phases`: the index m (folded to 0..N/2) of the m-twisted state psi_m that the ring is
close to, or -1 if it is not close to any. Uses Z = (1/N) sum_j exp(i (theta_{j+1} - theta_j)) (cyclic): the ring
counts as close to psi_m if |Z| > threshold, with arg Z = 2 pi m / N. psi_m and psi_{N-m} are mirror images.
"""
function twist_number(phases; threshold = 0.95)
    N = size(phases, 1)
    m = fill(-1, size(phases, 2))
    for k in axes(phases, 2)
        Z = sum(exp.(im .* (circshift(phases[:, k], -1) .- phases[:, k]))) / N
        if abs(Z) > threshold
            twist = mod(round(Int, angle(Z) * N / 2pi), N)
            m[k] = min(twist, N - twist)
        end
    end
    return m
end

"""
    write_omega_tau(prefix, times, phases, speeds)

Writes the projections of a trajectory onto the plane (tau = 1/v, frequency) used in Fig. 3(a),(b):
  <prefix>_omega_tau_mean.csv    columns mean_j(1/v_j), d/dt mean_j theta_j
  <prefix>_omega_tau_node1.csv   columns 1/v_1,          d/dt theta_1
(frequencies by forward differences between saved times).
"""
function write_omega_tau(prefix, times, phases, speeds)
    mean_phase = vec(mean(phases, dims = 1))
    mean_tau = [mean(1 ./ speeds[:, k]) for k in 2:length(times)]
    writedlm("$(prefix)_omega_tau_mean.csv", hcat(mean_tau, diff(mean_phase) ./ diff(times)), ',')
    writedlm("$(prefix)_omega_tau_node1.csv", hcat(1 ./ speeds[1, 2:end], diff(phases[1, :]) ./ diff(times)), ',')
    println("wrote $(prefix)_omega_tau_mean.csv and $(prefix)_omega_tau_node1.csv")
end
