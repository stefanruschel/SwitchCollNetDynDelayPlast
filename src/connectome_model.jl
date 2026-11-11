# connectome_model.jl
#
# The model of the Letter on the 68-node structural brain connectome (Figs. 1, 2 and S1):
#
#   theta_i'(t) = omega + sigma * amplitude * sum_j w_ij h( theta_j(t - tau_ij) - theta_i(t) ),      Eq. (1)
#   tau_ij      = d_ij / v_j,                                                                         Eq. (2)
#   v_j'(t)     = epsilon * g(theta_j', v_j)          (epsilon = 0 for t < t_on),
#
# with the interaction function
#
#   h(x) = -sin(x - a) + r sin(2x).
#
# Examples (the values of a are those used for the published runs):
#   a = pi,          r = 0               h(x) = sin x
#   a = pi + pi/4,   r = 0               h(x) = sin(x - pi/4)
#   a = pi + pi/2,   r = 1               h(x) = sin(x - pi/2) + sin 2x
# A phase lag alpha in the notation h(x) = sin(x - alpha) of the figure captions corresponds to a = alpha + pi.
# The prefactor `amplitude` (c_1, c_2 in Fig. 2) scales h so that max h = 1.
#
# Two adaptation rules g are implemented:
#
#   :sigmoid      Eq. (5):  g = -v_j + v_min + (v_max - v_min) f( F_j - nu ),   f(x) = 1 / (1 + exp(-beta x)),
#                 where F_j = (theta_j(t) - theta_j(t - window)) / window is the mean frequency of node j over
#                 the last `window` = 0.1 time units (memory kernel k = 10 * 1_[-0.1, 0] in the Letter).
#   :polynomial   Eq. (6):  g = (theta_j' - nu) (v_j - v_min) (v_max - v_j).
#
# Numerical method: DelayDiffEq.jl, MethodOfSteps(Tsit5()), with abstol = reltol = `tolerance`
# (chosen by each script; see the README for the values used for each figure).
#
# Implementation note: some arithmetic below is written in a particular order (for example the two-step
# normalisation of the weights). This is deliberate: it is the order used for the published runs, so that
# rerunning the scripts reproduces the published data files exactly rather than only up to rounding.

using DifferentialEquations, DelayDiffEq, LinearAlgebra, Random, Statistics, MAT

const CONNECTOME_DATA = joinpath(@__DIR__, "..", "connectome_data")
const NUMBER_OF_NODES = 68

# ------------------------------------------------------------------------------------------------ data

"""
    raw_connectome()

Returns the unnormalised matrices `(weights, distances)` of the 68-region Desikan-Killiany connectome:
structural weights `sc` from `sc_raw.mat` and Euclidean inter-regional distances `dist` from `D.mat`
(first slice of each file).
"""
function raw_connectome()
    weights = read_mat_variable("sc_raw.mat", "sc")
    distances = read_mat_variable("D.mat", "dist")
    return weights, distances
end

function read_mat_variable(filename, variable)
    file = matopen(joinpath(CONNECTOME_DATA, filename))
    matrix = read(file, variable)[:, :, 1]
    close(file)
    return matrix
end

"""
    load_connectome()

Returns `(W, D)` normalised as in the Letter: the weights rescaled to an average row sum of 1 (so that sigma is
the total input per node) and the distances rescaled to mean 1.
"""
function load_connectome()
    weights, distances = raw_connectome()
    # Average row sum 1, in two steps (first mean weight 1, then mean row sum 1) exactly as in the published runs.
    weights = weights ./ mean(weights)
    weights = weights ./ mean(sum(weights, dims = 2))
    distances = distances ./ mean(distances)
    return weights, distances
end

# ------------------------------------------------------------------------------------------------ model

"Interaction function h(x) = -sin(x - a) + r sin(2x)."
interaction_function(x, a, r) = -sin(x - a) + r * sin(2x)

"Sigmoid f(x) = 1 / (1 + exp(-beta x)) of the adaptation rule (5)."
sigmoid(x, beta) = 1 / (1 + exp(-beta * x))

"""
All parameters of the connectome model. Construct with keywords, e.g.

    ConnectomeParameters(W = W, D = D, sigma = 0.6, nu = 1.44, a = pi, r = 0.0, t_on = 1000.0)

Fields without a default must be given.
"""
Base.@kwdef struct ConnectomeParameters
    W::Matrix{Float64}            # coupling weights w_ij
    D::Matrix{Float64}            # distances d_ij
    sigma::Float64                # coupling strength
    nu::Float64                   # threshold of the adaptation rule
    a::Float64                    # phase parameter of h
    r::Float64                    # weight of the second harmonic of h
    amplitude::Float64 = 1.0      # prefactor of h (c_1, c_2 in Fig. 2)
    omega::Float64 = 1.0          # natural frequency (identical for all nodes)
    rule::Symbol = :sigmoid       # adaptation rule, :sigmoid (Eq. (5)) or :polynomial (Eq. (6))
    epsilon::Float64 = 0.01       # adaptation rate
    t_on::Float64 = 0.0           # adaptation is switched on at t = t_on (speeds constant before)
    beta::Float64 = 5.0           # gain of the sigmoid
    window::Float64 = 0.1         # averaging window of the sigmoid rule (memory kernel 10 * 1_[-0.1, 0])
    vmin::Float64 = 0.1           # lower bound of the speeds
    vmax::Float64 = 20.0          # upper bound of the speeds
end

"Largest delay that can occur (distance / smallest speed); needed by the solver to store enough history."
maximum_delay(p::ConnectomeParameters) = maximum(p.D) / p.vmin

"""
    connectome_rhs!(du, u, history, p, t)

Right-hand side of the delay differential equation. The state is u = (theta_1..theta_N, v_1..v_N).
`history(p, s; idxs = j)` returns theta_j(s) at an earlier time s.
"""
function connectome_rhs!(du, u, history, p::ConnectomeParameters, t)
    N = size(p.W, 1)
    theta = view(u, 1:N)
    speed = view(u, N+1:2N)
    dtheta = view(du, 1:N)
    dspeed = view(du, N+1:2N)

    # phases, Eq. (1)
    for i in 1:N
        coupling = 0.0
        for j in 1:N
            if p.W[i, j] != 0.0
                delay = p.D[i, j] / max(speed[j], p.vmin)
                delayed_phase = history(p, t - delay; idxs = j)
                coupling += p.W[i, j] * interaction_function(delayed_phase - theta[i], p.a, p.r)
            end
        end
        dtheta[i] = p.omega + p.sigma * p.amplitude * coupling
    end

    # speeds, Eq. (5) or (6)
    epsilon = t < p.t_on ? 0.0 : p.epsilon
    if p.rule == :sigmoid
        for j in 1:N
            earlier_phase = history(p, t - p.window; idxs = j)
            drive = (theta[j] - earlier_phase) / p.window - p.nu      # mean frequency over [t - window, t] minus threshold
            dspeed[j] = epsilon * (-speed[j] + p.vmin + (p.vmax - p.vmin) * sigmoid(drive, p.beta))
        end
    elseif p.rule == :polynomial
        for j in 1:N
            dspeed[j] = epsilon * (dtheta[j] - p.nu) * (speed[j] - p.vmin) * (p.vmax - speed[j])
        end
    else
        error("unknown adaptation rule $(p.rule)")
    end
    return nothing
end

"""
    initial_history(theta0, v0, omega)

History for t <= 0: the phases rotate with the natural frequency from theta0, theta_j(t) = theta0_j + omega t,
and the speeds are constant, v_j(t) = v0_j. Supports the `idxs` keyword used by the solver.
"""
function initial_history(theta0, v0, omega)
    N = length(theta0)
    function history(p, t; idxs = nothing)
        if idxs === nothing
            return vcat(theta0 .+ omega .* t, v0)
        elseif idxs isa Integer
            return idxs <= N ? theta0[idxs] + omega * t : v0[idxs-N]
        else
            return [k <= N ? theta0[k] + omega * t : v0[k-N] for k in idxs]
        end
    end
    return history
end

"""
    random_initial_condition(rng, N; vmin, vmax, mean_speed = 2.5, speed_spread = 0.1)

theta_i ~ U(0, 2 pi) and v_i ~ N(mean_speed, speed_spread) restricted to [vmin, vmax], drawn from `rng`.
"""
function random_initial_condition(rng, N; vmin, vmax, mean_speed = 2.5, speed_spread = 0.1)
    theta0 = 2pi .* rand(rng, N)
    v0 = clamp.(mean_speed .+ speed_spread .* randn(rng, N), vmin, vmax)
    return theta0, v0
end

"""
    simulate(p, theta0, v0; tend, saveat, tolerance)

Integrates the model from t = 0 to `tend` and returns the solution object (saved at the times `saveat`).
"""
function simulate(p::ConnectomeParameters, theta0, v0; tend, saveat, tolerance)
    history = initial_history(theta0, v0, p.omega)
    u0 = vcat(theta0, v0)
    problem = DDEProblem(connectome_rhs!, u0, history, (0.0, tend), p; constant_lags = [p.window, maximum_delay(p)])
    return solve(problem, MethodOfSteps(Tsit5()); abstol = tolerance, reltol = tolerance, saveat = saveat)
end

# ------------------------------------------------------------------------------------------------ diagnostics

"""
    solution_window(solution, N, tmin, tmax)

Saved times in [tmin, tmax] and the corresponding phases (N x n matrix) and speeds (N x n matrix).
"""
function solution_window(solution, N, tmin, tmax)
    indices = findall(t -> tmin <= t <= tmax, solution.t)
    states = reduce(hcat, solution.u[indices])
    return solution.t[indices], states[1:N, :], states[N+1:2N, :]
end

"""
    order_parameter(phases, m = 1)

R_m(t) = |sum_j exp(i m theta_j(t))| / N for every column (time) of the N x n matrix `phases`.
R_1 is the Kuramoto order parameter, R_2 measures two-cluster states.
"""
order_parameter(phases, m = 1) = [abs(sum(exp.(im * m .* phases[:, k]))) / size(phases, 1) for k in 1:size(phases, 2)]

"Mean frequency of every node over the time span of `times` (from the unwrapped phases)."
mean_frequencies(times, phases) = (phases[:, end] .- phases[:, 1]) ./ (times[end] - times[1])

"Fraction of nodes whose mean frequency lies within `tolerance` of the median mean frequency."
fraction_locked(frequencies; tolerance = 1e-3) = count(abs.(frequencies .- median(frequencies)) .< tolerance) / length(frequencies)

"""
    count_collapses(R; high, low)

Number of collapses of the order parameter time series R: R rises above `high` and later falls below `low`
(counted once per rise). With high = 0.9, low = 0.1 this is the switching criterion of Fig. 2.
"""
function count_collapses(R; high, low)
    collapses = 0
    armed = false
    for value in R
        if value > high
            armed = true
        end
        if armed && value < low
            collapses += 1
            armed = false
        end
    end
    return collapses
end

"""
    bell_order(frequencies)

Permutation of the nodes that puts the fastest node in the middle and slower nodes alternately below and above it
(used to reorder the raster plot of Fig. 1(b2) and the frequencies of Fig. 1(b3)).
"""
function bell_order(frequencies)
    N = length(frequencies)
    descending = sortperm(frequencies, rev = true)
    order = zeros(Int, N)
    middle = (N + 1) ÷ 2
    order[middle] = descending[1]
    k = 2
    for d in 1:(middle-1)
        if k <= N
            order[middle-d] = descending[k]
            k += 1
        end
        if k <= N
            order[middle+d] = descending[k]
            k += 1
        end
    end
    if iseven(N)
        order[N] = descending[N]
    end
    return order
end
