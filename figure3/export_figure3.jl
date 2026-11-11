# export_figure3.jl
#
# Writes the plain-text tables that Figure3.tex (pgfplots) draws. Unlike the other export scripts this one also
# integrates the blue trajectory of Fig. 3 ("switching dynamics", about one minute), because panels (a), (b) and (c)
# are all drawn from that run.
#
# Usage (from this folder):
#   julia --project=.. export_figure3.jl        then        lualatex Figure3.tex
#
# Inputs:
#   data/twisted_state_branches.mat                branches Omega_m(tau) of the m-twisted states psi_m of the ring at
#                                                  epsilon = 0 and the number of unstable eigenvalues along them
#                                                  (twisted_state_branches.m, MATLAB + DDE-Biftool)
#   data/frequency_selection_omega_tau_node1.csv   green trajectory (simulate_frequency_selection.jl)
# Outputs (tables/):
#   branch_m<m>_stable.dat, branch_m<m>_unstable.dat   branches, split into stable and unstable parts
#   bif_m<m>.dat                                       points where the stability changes (markers)
#   eq_nu.dat, eq_bounds.dat                           equilibria of the slow flow: branch points with Omega = nu, and
#                                                      branch points on the bounds tau = 1/v_max and tau = 1/v_min
#   traj_blue_node1.dat, traj_blue_mean.dat            blue trajectory projected onto (1/v, frequency): node 1 and mean
#   traj_green.dat                                     green trajectory (node 1)
#   jumps_b.dat                                        arrows for the fast jumps in panel (b)
#   panel_c.dat, labels_c.dat                          R(t), C(t) and the twisted state of each segment for panel (c)
#   params.tex                                         time range of panel (c)
#
# The blue trajectory: theta_j(t) = Omega t, v_j(t) = c for t <= 0 with Omega = 1.0004 and c = 1/5 (on the synchronous
# branch psi_0, just above the threshold nu = 1), integrated to t = 2.4e6 so that one full cycle psi_0 -> ... -> psi_0
# is contained. Panel (c) shows this cycle, from 0.06e6 before the departure from psi_0 to 0.06e6 after the return.

include(joinpath(@__DIR__, "..", "src", "ring_model.jl"))
using MAT, DelimitedFiles, Printf

const TABLES = joinpath(@__DIR__, "tables")
const N = 11
const NU = 1.0
const TAU_MIN = 1 / 5.0       # 1/v_max
const TAU_MAX = 6.8           # 1/v_min

"Writes rows of numbers separated by spaces; a row starting with NaN becomes 'nan nan' (break between segments)."
function write_rows(name, rows; header = "x y")
    open(joinpath(TABLES, name), "w") do io
        println(io, header)
        for row in rows
            if isnan(row[1])
                println(io, "nan nan")
            else
                println(io, join((@sprintf("%.6g", value) for value in row), ' '))
            end
        end
    end
end

"""
Branches of the twisted states from the MATLAB results: stable and unstable parts, stability changes, and the
equilibria of the slow flow (Omega = nu within the bounds, and branch points on the bounds).
"""
function export_branches()
    file = matopen(joinpath(@__DIR__, "data", "twisted_state_branches.mat"))
    Omega_all = read(file, "Omega")
    tau_all = read(file, "tau1")
    unstable_all = read(file, "Npos")       # number of unstable eigenvalues
    close(file)

    equilibria_nu = Vector{Float64}[]
    equilibria_bounds = Vector{Float64}[]
    for matlab_index in (11, 1, 2, 3, 4, 5)   # MATLAB stores psi_0 under index N = 11
        m = matlab_index % N
        Omega = vec(Omega_all[matlab_index])
        tau = vec(tau_all[matlab_index])
        stable = vec(unstable_all[matlab_index]) .== 0

        for (name, wanted) in (("stable", true), ("unstable", false))
            rows = Vector{Float64}[]
            for i in eachindex(tau)
                if stable[i] == wanted
                    push!(rows, [tau[i], Omega[i]])
                elseif i > 1 && stable[i-1] == wanted       # close the segment at the change of stability, then break
                    push!(rows, [tau[i], Omega[i]])
                    push!(rows, [NaN, NaN])
                end
            end
            write_rows("branch_m$(m)_$(name).dat", rows)
        end

        changes = [[tau[i], Omega[i]] for i in 2:length(tau) if stable[i] != stable[i-1] && tau[i] > 0]
        write_rows("bif_m$(m).dat", changes)

        for i in 1:length(tau)-1
            # crossing of Omega = nu between the bounds
            if (Omega[i] - NU) * (Omega[i+1] - NU) < 0
                s = (NU - Omega[i]) / (Omega[i+1] - Omega[i])
                t = tau[i] + s * (tau[i+1] - tau[i])
                if TAU_MIN <= t <= TAU_MAX
                    push!(equilibria_nu, [t, NU])
                end
            end
            # crossing of the bound lines tau = 1/v_max and tau = 1/v_min
            for bound in (TAU_MIN, TAU_MAX)
                if (tau[i] - bound) * (tau[i+1] - bound) <= 0 && tau[i] != tau[i+1]
                    s = (bound - tau[i]) / (tau[i+1] - tau[i])
                    push!(equilibria_bounds, [bound, Omega[i] + s * (Omega[i+1] - Omega[i])])
                end
            end
        end
    end
    write_rows("eq_nu.dat", equilibria_nu)
    write_rows("eq_bounds.dat", equilibria_bounds)
end

"""
Visits of twisted states: runs of constant twist number m >= 0 lasting longer than `minimum_length`,
as rows [start time, end time, m].
"""
function visits(times, twist; minimum_length = 1.5e4)
    result = Vector{Float64}[]
    k = 1
    while k <= length(times)
        if twist[k] >= 0
            j = k
            while j < length(times) && twist[j+1] == twist[k]
                j += 1
            end
            if times[j] - times[k] > minimum_length
                push!(result, [times[k], times[j], twist[k]])
            end
            k = j + 1
        else
            k += 1
        end
    end
    return result
end

function main()
    mkpath(TABLES)
    export_branches()

    # green trajectory (frequency selection): node 1
    green = readdlm(joinpath(@__DIR__, "data", "frequency_selection_omega_tau_node1.csv"), ',')
    write_rows("traj_green.dat", [green[i, :] for i in 1:size(green, 1)])

    # blue trajectory (switching dynamics)
    s = RingSystem(RingParameters())
    tend = 2.4e6
    save_times = vcat(collect(0.0:1000.0:1.0e6), collect(1.0e6+10:10.0:tend))
    elapsed = @elapsed solution = simulate_from_locked_state(s, 1.0004, 1 / 5; tend = tend, saveat = save_times)
    println("blue trajectory integrated to t = $tend in $(round(elapsed, digits = 1)) s")
    times, phases, speeds = solution_window(solution, N, 0.0, tend)

    # one cycle: from the end of the initial visit of psi_0 to the start of the next visit of psi_0
    state_visits = visits(times, twist_number(phases))
    departure = state_visits[1][2]
    next_psi0 = findfirst(i -> state_visits[i][3] == 0, 2:length(state_visits))
    if next_psi0 === nothing
        error("no return to psi_0 before t = $tend")
    end
    return_time = state_visits[next_psi0+1][1]
    t0 = departure - 0.06e6
    t1 = min(return_time + 0.06e6, tend)
    @printf("departure from psi_0 at t = %.4g, return at t = %.4g; panel (c) shows [%.4g, %.4g]\n", departure, return_time, t0, t1)
    println("visits: ", join([@sprintf("psi_%d", Int(v[3])) for v in state_visits if v[2] >= t0 && v[1] <= t1], " -> "))

    # (tau, frequency) projections up to the end of panel (c), every 4th saved time
    mean_phase = vec(mean(phases, dims = 1))
    selected = findall(<=(t1), times)
    rows_mean = [[mean(1 ./ speeds[:, k]), (mean_phase[k] - mean_phase[k-1]) / (times[k] - times[k-1])] for k in selected[2:4:end]]
    rows_node1 = [[1 / speeds[1, k], (phases[1, k] - phases[1, k-1]) / (times[k] - times[k-1])] for k in selected[2:4:end]]
    write_rows("traj_blue_mean.dat", rows_mean)
    write_rows("traj_blue_node1.dat", rows_node1)

    # fast jumps of the mean frequency across Omega = 0.96 (between the psi_5 branch and the upper branches) for the
    # blue arrows of panel (b): first crossing at each tau position (jumps closer than 0.02 in tau are merged)
    jumps = Vector{Float64}[]
    for k in 2:length(rows_mean)
        before, after = rows_mean[k-1][2], rows_mean[k][2]
        if (before - 0.96) * (after - 0.96) < 0 && rows_mean[k][1] > 2.9
            direction = after < before ? -1.0 : 1.0
            if !any(j -> abs(j[1] - rows_mean[k][1]) < 0.02, jumps)
                push!(jumps, [rows_mean[k][1], direction])
            end
        end
    end
    # arrows beside the jump (left of downward, right of upward jumps): position (x, y) and direction (u, v)
    arrows = [[j[1] - 0.012 * j[2], j[2] < 0 ? 1.0 : 0.94, 0.0, j[2] * 0.05] for j in jumps]
    write_rows("jumps_b.dat", arrows; header = "x y u v")
    println("fast jumps (tau, direction): ", join([@sprintf("(%.3f,%+d)", j[1], Int(j[2])) for j in jumps], " "))

    # panel (c): order parameter R and nearest-neighbour coherence C (every 3rd saved time), time in units of 1e6
    window = findall(t -> t0 <= t <= t1, times)[1:3:end]
    R = order_parameter(phases[:, window])
    C = nearest_neighbour_coherence(phases[:, window])
    write_rows("panel_c.dat", [[times[window[k]] / 1e6, R[k], C[k]] for k in eachindex(window)]; header = "t R C")
    # label psi_m at the centre of the visible part of every visit longer than 0.02e6, at alternating heights
    labels = [[(max(v[1], t0) + min(v[2], t1)) / 2e6, v[3]] for v in state_visits if min(v[2], t1) - max(v[1], t0) > 0.02e6]
    write_rows("labels_c.dat", [[l[1], l[2], isodd(i) ? 0.62 : 0.38] for (i, l) in enumerate(labels)]; header = "t m y")

    open(joinpath(TABLES, "params.tex"), "w") do io
        @printf(io, "\\def\\CtMin{%.4f}\\def\\CtMax{%.4f}\n", t0 / 1e6, t1 / 1e6)
    end
    println("wrote tables/*")
end

main()
