# run_figure2_sweeps.jl
#
# Produces the data of Fig. 2 of the Letter and the neighbourhood checks of the Fig. 1 parameters (SM, "Sweeps of Fig. 2").
#
# Usage (from this folder; -t auto uses all cores):
#   julia -t auto --project=.. run_figure2_sweeps.jl                  the six 45 x 45 sweeps of Fig. 2(b)-(g)
#   julia -t auto --project=.. run_figure2_sweeps.jl --only=b,g       only some panels
#   julia -t auto --project=.. run_figure2_sweeps.jl neighbourhoods   3 x 3 cells x 4 seeds around the Fig. 1 points
#   julia -t auto --project=.. run_figure2_sweeps.jl --quick          3 x 3 grid, t <= 200 (test, files quick_*.csv)
#
# Run time: about 2025 runs of 5-25 s each per panel, i.e. roughly 3 hours for all six panels on 14 threads.
#
# Outputs (data/):
#   sweep_figure2<panel>.csv          panel = b..g, see INTERACTION_FUNCTIONS in sweep_connectome.jl
#   neighbourhood_figure1b.csv        chimera of Fig. 1(b), h = sin x                (around sigma = 0.6, nu = 1.44)
#   neighbourhood_figure1d.csv        cluster switching of Fig. 1(d), panel (g)      (around sigma = 1.0, nu = 1.4)
#   neighbourhood_figure2f.csv        a switching point in panel (f), not shown in Fig. 1
# Columns: see sweep_connectome.jl.

include(joinpath(@__DIR__, "sweep_connectome.jl"))
include(joinpath(@__DIR__, "..", "src", "command_line.jl"))

# Grid of Fig. 2: sigma in [0.3, 3] (45 values, logarithmic), nu in [1.0, 1.6] (45 values), one initial condition
# (seed 1) per cell, t <= 1500.
const FIGURE2_GRID = (sigma_min = 0.3, sigma_max = 3.0, n_sigma = 45, nu_min = 1.0, nu_max = 1.6, n_nu = 45,
                      n_seeds = 1, tend = 1500.0)

# Neighbourhood checks: 3 x 3 cells around a point, 4 initial conditions each, t <= 2000.
const NEIGHBOURHOODS = [
    (file = "neighbourhood_figure1b.csv", panel = "b", sigma_min = 0.5, sigma_max = 0.72, nu_min = 1.40, nu_max = 1.48),
    (file = "neighbourhood_figure1d.csv", panel = "g", sigma_min = 0.9, sigma_max = 1.17, nu_min = 1.37, nu_max = 1.45),
    (file = "neighbourhood_figure2f.csv", panel = "f", sigma_min = 1.7, sigma_max = 2.05, nu_min = 1.23, nu_max = 1.31),
]

function main()
    data = joinpath(@__DIR__, "data")
    mkpath(data)

    if "neighbourhoods" in ARGS
        for check in NEIGHBOURHOODS
            run_sweep(joinpath(data, check.file); interaction = INTERACTION_FUNCTIONS[check.panel],
                      sigma_min = check.sigma_min, sigma_max = check.sigma_max, n_sigma = 3,
                      nu_min = check.nu_min, nu_max = check.nu_max, n_nu = 3, n_seeds = 4, tend = 2000.0)
        end
        return
    end

    panels = split(option("only", "b,c,d,e,f,g"), ",")
    quick = flag("quick")
    grid = quick ? merge(FIGURE2_GRID, (n_sigma = 3, n_nu = 3, tend = 200.0)) : FIGURE2_GRID
    for panel in panels
        file = joinpath(data, (quick ? "quick_" : "") * "sweep_figure2$(panel).csv")
        run_sweep(file; interaction = INTERACTION_FUNCTIONS[panel],
                  sigma_min = grid.sigma_min, sigma_max = grid.sigma_max, n_sigma = grid.n_sigma,
                  nu_min = grid.nu_min, nu_max = grid.nu_max, n_nu = grid.n_nu, n_seeds = grid.n_seeds, tend = grid.tend)
    end
end

main()
