# export_figureS1.jl
#
# Writes the plain-text tables that FigureS1.tex (pgfplots) draws. No simulation.
#
# Usage (from this folder):
#   julia --project=.. export_figureS1.jl        then        lualatex FigureS1.tex
#
# Inputs:  data/figureS1a.csv      (simulate_figureS1a.jl)
#          data/sweep1d_*.csv      (run_one_parameter_sweeps.jl)
# Outputs: tables/speeds_a.dat     speeds every 10 time units          (panel (a))
#          tables/order_a.dat      R1 and R2 every 2 time units        (panel (a))
#          tables/share_*.dat      share of switching runs per parameter value, columns x (value) and s (share)
#                                  (panels (b)-(f)); a run switches if R1 or R2 collapses (0.9 -> 0.1) at least twice

include(joinpath(@__DIR__, "..", "src", "connectome_model.jl"))
using DelimitedFiles, Printf

const TABLES = joinpath(@__DIR__, "tables")
data_file(name) = joinpath(@__DIR__, "data", name)

"Writes a table with a header line, numbers rounded to 5 significant digits."
function write_table(name, header, matrix)
    open(joinpath(TABLES, name), "w") do io
        println(io, header)
        writedlm(io, round.(matrix; sigdigits = 5), ' ')
    end
end

"""
    switching_share(files...)

For every parameter value in the sweep files: the fraction of runs that switch. If a value occurs in several files,
the first file wins (so a refined sweep can be listed before the coarse one). Returns a matrix [value share], sorted.
"""
function switching_share(files...)
    values = Float64[]
    shares = Float64[]
    for file in files
        table, header = readdlm(data_file(file), ','; header = true)
        header = vec(header)
        column(name) = table[:, findfirst(==(name), header)]
        for value in unique(column("value"))
            runs = column("value") .== value
            switches = (column("events_R1s")[runs] .>= 2) .| (column("events_R2s")[runs] .>= 2)
            push!(values, value)
            push!(shares, mean(switches))
        end
    end
    order = sortperm(values)
    first_occurrences = unique(i -> values[i], order)
    return hcat(values[first_occurrences], shares[first_occurrences])
end

function main()
    mkpath(TABLES)

    # (a) polynomial rule (6): speeds and order parameters
    table, _ = readdlm(data_file("figureS1a.csv"), ','; header = true)
    N = NUMBER_OF_NODES
    times = table[:, 1]
    phases = permutedims(table[:, 2:N+1])
    speeds = permutedims(table[:, N+2:end])
    speed_indices = 1:10:length(times)
    order_indices = 1:2:length(times)
    write_table("speeds_a.dat", join(vcat("t", ["v$i" for i in 1:N]), ' '), hcat(times[speed_indices], permutedims(speeds[:, speed_indices])))
    write_table("order_a.dat", "t R1 R2", hcat(times[order_indices], order_parameter(phases[:, order_indices], 1),
                                               order_parameter(phases[:, order_indices], 2)))

    # (b)-(f) share of switching runs
    write_table("share_eps.dat", "x s", switching_share("sweep1d_eps.csv"))
    write_table("share_spread.dat", "x s", switching_share("sweep1d_spread.csv"))
    write_table("share_beta.dat", "x s", switching_share("sweep1d_beta.csv"))
    write_table("share_nu_beta3.dat", "x s", switching_share("sweep1d_nu_beta3.csv"))
    write_table("share_nu_beta8.dat", "x s", switching_share("sweep1d_nu_beta8.csv"))
    write_table("share_rule6_eps0.01.dat", "x s", switching_share("sweep1d_rule6_nu_fine.csv", "sweep1d_rule6_nu.csv"))
    write_table("share_rule6_eps0.002.dat", "x s", switching_share("sweep1d_rule6_nu_eps0.002_fine.csv", "sweep1d_rule6_nu_eps0.002.csv"))

    for file in sort(readdir(TABLES))
        if startswith(file, "share")
            shares = readdlm(joinpath(TABLES, file), ' '; skipstart = 1)[:, 2]
            println(file, ": ", join(string.(round.(shares; digits = 2)), " "))
        end
    end
    println("wrote tables/*")
end

main()
