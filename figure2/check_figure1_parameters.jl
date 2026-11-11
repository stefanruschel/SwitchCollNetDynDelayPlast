# check_figure1_parameters.jl
#
# Checks, from the sweep data, the statements in the SM (paragraph "Sweeps of Fig. 2") on how robust the parameters
# of the Fig. 1 time series are. Prints the number of qualifying runs in every cell of the 3 x 3 neighbourhoods
# (4 initial conditions per cell) and the classification of the 5 x 5 block of the Fig. 2(b) grid:
#   Fig. 1(b)  chimera-like at the chosen point (sigma, nu) = (0.6, 1.44), the centre of its neighbourhood
#   Fig. 1(d)  switching at the centre of its neighbourhood, (1.026, 1.41), next to the chosen point (1.0, 1.4)
#   Fig. 1(c)  all 24 cells of the 5 x 5 block of the 45 x 45 grid of Fig. 2(b) around (0.8, 1.56) show switching
#
# Usage (from this folder):  julia --project=.. check_figure1_parameters.jl
# Inputs: data/neighbourhood_figure1b.csv, data/neighbourhood_figure1d.csv, data/sweep_figure2b.csv

using DelimitedFiles, Statistics, Printf

is_switching(collapses_R1, collapses_R2) = collapses_R1 >= 2 || collapses_R2 >= 2
is_chimera_like(fraction_locked, R1_mean) = 0.1 < fraction_locked < 0.9 && 0.2 < R1_mean < 0.9

function read_sweep(file)
    table, header = readdlm(joinpath(@__DIR__, "data", file), ','; header = true)
    header = vec(header)
    column(name) = table[:, findfirst(==(name), header)]
    return (sigma = column("sigma"), nu = column("nu"), R1 = column("R1mean"), locked = column("frac_locked"),
            collapses_R1 = column("events_R1_strict"), collapses_R2 = column("events_R2_strict"))
end

"Prints, for every cell of a neighbourhood sweep, how many runs satisfy `criterion`."
function report_neighbourhood(file, description, criterion)
    s = read_sweep(file)
    println("$description ($file):")
    all_runs = 0
    passing_runs = 0
    for sigma in sort(unique(s.sigma)), nu in sort(unique(s.nu))
        runs = findall((s.sigma .== sigma) .& (s.nu .== nu))
        passing = count(k -> criterion(s, k), runs)
        @printf("  sigma = %.3f  nu = %.3f   %d of %d runs\n", sigma, nu, passing, length(runs))
        all_runs += length(runs)
        passing_runs += passing
    end
    println("  total: $passing_runs of $all_runs runs\n")
end

report_neighbourhood("neighbourhood_figure1b.csv", "Fig. 1(b), chimera-like", (s, k) -> is_chimera_like(s.locked[k], s.R1[k]))
report_neighbourhood("neighbourhood_figure1d.csv", "Fig. 1(d), switching", (s, k) -> is_switching(s.collapses_R1[k], s.collapses_R2[k]))

# Fig. 1(c): 5 x 5 block of the Fig. 2(b) grid around the cell closest to (sigma, nu) = (0.8, 1.56)
s = read_sweep("sweep_figure2b.csv")
sigmas = sort(unique(s.sigma))
nus = sort(unique(s.nu))
j0 = argmin(abs.(log.(sigmas) .- log(0.8)))
i0 = argmin(abs.(nus .- 1.56))
@printf("Fig. 1(c), switching in the 5 x 5 block of Fig. 2(b) around the cell (sigma, nu) = (%.3f, %.3f):\n", sigmas[j0], nus[i0])
neighbours = 0
switching_neighbours = 0
for i in i0-2:i0+2, j in j0-2:j0+2
    (i == i0 && j == j0) && continue
    k = findfirst((s.sigma .== sigmas[j]) .& (s.nu .== nus[i]))
    global neighbours += 1
    global switching_neighbours += is_switching(s.collapses_R1[k], s.collapses_R2[k])
end
k0 = findfirst((s.sigma .== sigmas[j0]) .& (s.nu .== nus[i0]))
println("  centre cell switching: $(is_switching(s.collapses_R1[k0], s.collapses_R2[k0]))")
println("  surrounding cells switching: $switching_neighbours of $neighbours")
