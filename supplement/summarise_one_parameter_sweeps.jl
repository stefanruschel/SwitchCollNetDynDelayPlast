# summarise_one_parameter_sweeps.jl
#
# Prints, for every parameter value of the one-parameter sweeps, how many runs switch and some averages; these are
# the numbers quoted in the SM paragraph "One-parameter sweeps".
#
# Usage (from this folder):
#   julia --project=.. summarise_one_parameter_sweeps.jl                         all data/sweep1d_*.csv
#   julia --project=.. summarise_one_parameter_sweeps.jl data/sweep1d_eps.csv    selected files
#
# Columns of the output:
#   value      value of the swept parameter
#   switching  number of runs with at least two strict collapses (0.9 -> 0.1) of R1 or R2, out of all runs
#   <R2 coll.> mean number of strict collapses of R2
#   <R1>, <R2> mean time-averaged order parameters
#   period     (run length / 2) / <R2 coll.>, a rough switching period (statistics are over the second half of the run)

using DelimitedFiles, Statistics, Printf

function summarise(file)
    table, header = readdlm(file, ','; header = true)
    header = vec(header)
    column(name) = table[:, findfirst(==(name), header)]
    println("--- $file")
    println("    value  switching  <R2 coll.>   <R1>   <R2>   period")
    for value in sort(unique(column("value")))
        runs = column("value") .== value
        switching = (column("events_R2s")[runs] .>= 2) .| (column("events_R1s")[runs] .>= 2)
        run_length = column("tend")[runs][1]
        collapses = mean(column("events_R2s")[runs])
        period = collapses > 0 ? @sprintf("%.0f", run_length / 2 / collapses) : "-"
        @printf("%9.4g   %d of %d    %8.2f   %5.2f  %5.2f   %s\n", value, count(switching), count(runs), collapses,
                mean(column("R1mean")[runs]), mean(column("R2mean")[runs]), period)
    end
end

files = isempty(ARGS) ? sort(filter(f -> startswith(basename(f), "sweep1d_"), readdir(joinpath(@__DIR__, "data"); join = true))) : ARGS
foreach(summarise, files)
