# Switching between Collective Network Dynamics driven by Time-Delay Plasticity

Code and data for

> S. Ruschel, E. Hristov, H. G. E. Meijer, S. Coombes, and R. Nicks,
> *Switching between Collective Network Dynamics driven by Time-Delay Plasticity*, [arXiv:2605.23520](https://arxiv.org/abs/2605.23520).

Repository: <https://github.com/stefanruschel/SwitchCollNetDynDelayPlast>

The scripts reproduce every numerical result of the Letter and its Supplemental Material (SM): the simulations, the
parameter sweeps, the tables of the SM and the figures themselves. All simulation results are included, so the
figures can be redrawn in a few minutes without rerunning any simulation.

## Contents

| Folder | Content |
|---|---|
| `src/` | shared code: the model on the connectome (`connectome_model.jl`), the model on the ring (`ring_model.jl`), command-line helpers |
| `connectome_data/` | the 68-region structural connectome (weights and distances), see `connectome_data/README.md` |
| `figure1/` | Fig. 1: chimera state, switching between synchrony and incoherence, switching between two-cluster states |
| `figure2/` | Fig. 2: sweeps over coupling strength and threshold for six interaction functions; robustness of the Fig. 1 parameters |
| `figure3/` | Fig. 3: ring of 11 nodes, branches of twisted states (MATLAB + DDE-Biftool) and the two trajectories |
| `supplement/` | SM: Fig. S1 (one-parameter sweeps, polynomial rule), Table S1 (ring, dependence on epsilon), persistence of the chimera state |

Every figure folder has the same structure:

| | |
|---|---|
| `simulate_*.jl`, `run_*.jl` | simulations; write their results to `data/` |
| `data/` | simulation results (CSV, included) |
| `export_*.jl` | read `data/` and write the plain-text tables that the figure draws to `tables/` (no simulation) |
| `tables/` | these tables (included) |
| `Figure*.tex` | the figure, drawn with pgfplots at its final size in the Letter |
| `Figure*.pdf` | the compiled figure |

Every script starts with a comment block that explains what it computes, how to call it, how long it takes and what
the columns of its output files mean.

## Requirements

* [Julia](https://julialang.org) 1.13 (tested with 1.13.1). The packages and their exact versions are listed in
  `Project.toml` and `Manifest.toml` (created with Julia 1.13.1; with another Julia version run `Pkg.resolve()` after
  activating the environment); install them once with

  ```
  julia setup.jl
  ```

* A LaTeX distribution with `lualatex`, `pgfplots` (1.18) and `standalone`, to compile the figures.
* Only for the branches of twisted states in Fig. 3: MATLAB and [DDE-Biftool](https://github.com/DDE-BifTool/DDE-Biftool)
  3.1. Their result is included (`figure3/data/twisted_state_branches.mat`), so MATLAB is not needed to redraw Fig. 3.

All commands below are run from the folder of the respective figure. `--project=..` selects the Julia environment of
this repository, and `-t auto` lets the sweeps use all processor cores.

## Redrawing the figures from the included data (a few minutes)

```
cd figure1;    julia --project=.. export_figure1.jl;   lualatex Figure1.tex
cd figure2;    julia --project=.. export_figure2.jl;   lualatex Figure2.tex
cd figure3;    julia --project=.. export_figure3.jl;   lualatex Figure3.tex     # also integrates the blue trajectory (~1 min)
cd supplement; julia --project=.. export_figureS1.jl;  lualatex FigureS1.tex
```

## Reproducing the simulations

Run times were measured on a Windows 11 desktop with 14 threads; the first call of each script additionally spends a
few minutes compiling.

### Fig. 1 (`figure1/`)

```
julia --project=.. simulate_figure1.jl b      # row (b), chimera state                                   ~1 min
julia --project=.. simulate_figure1.jl c      # row (c), switching between synchrony and incoherence    ~1 min
julia --project=.. simulate_figure1.jl d      # row (d), switching between two-cluster states            ~1 min
```

Output: `data/figure1<row>.csv` (time, unwrapped phases and speeds of all 68 nodes every 0.5 time units up to
t = 2000) and `data/figure1<row>.params.txt`.

### Fig. 2 (`figure2/`)

```
julia -t auto --project=.. run_figure2_sweeps.jl                   # six 45 x 45 sweeps, panels (b)-(g)    ~3 h
julia -t auto --project=.. run_figure2_sweeps.jl neighbourhoods    # checks of the Fig. 1 parameters       ~5 min
julia --project=.. check_figure1_parameters.jl                     # prints the robustness statements of the SM
```

`run_figure2_sweeps.jl --only=b,g` runs selected panels, `--quick` a 3 x 3 test grid. Output:
`data/sweep_figure2<panel>.csv` with one line per run (time-averaged order parameters, number of collapses of R1 and
R2, fraction of frequency-locked nodes, ...; the columns are explained in `sweep_connectome.jl`), and
`data/neighbourhood_*.csv` in the same format.

### Fig. 3 (`figure3/`)

```
julia --project=.. simulate_frequency_selection.jl     # green trajectory, t <= 3e6             ~1 min
julia --project=.. export_figure3.jl                   # blue trajectory (t <= 2.4e6) and tables  ~1 min
matlab -batch twisted_state_branches                   # branches and their stability            ~1 h
```

For the MATLAB script set the path to DDE-Biftool in its first lines.

### Supplemental Material (`supplement/`)

```
julia --project=.. simulate_figureS1a.jl                   # Fig. S1(a), polynomial rule (6), t <= 1e4   ~1 min
julia -t auto --project=.. run_one_parameter_sweeps.jl     # Fig. S1(b)-(f), nine sweeps                 ~1 h
julia --project=.. summarise_one_parameter_sweeps.jl       # numbers quoted in the SM text
julia -t auto --project=.. ring_epsilon_sweep.jl           # Table S1                                    ~1 min
julia --project=.. chimera_long_run.jl                     # persistence of the chimera, t <= 1e6        ~80 min
```

## The model and its parameters

The phases theta_i and the axonal speeds v_j of N nodes obey (Letter, Eqs. (1), (2), (5), (6))

```
theta_i'(t) = omega + sigma * sum_j w_ij h( theta_j(t - tau_ij) - theta_i(t) ),     tau_ij = d_ij / v_j,
v_j'(t)     = epsilon g(theta_j', v_j),
```

with the sigmoid rule (5), `g = -v_j + v_min + (v_max - v_min) / (1 + exp(-beta (F_j - nu)))`, where F_j is the mean
frequency of node j over the last 0.1 time units (memory kernel `k = 10 * 1_[-0.1, 0]`), or the polynomial rule (6),
`g = (theta_j' - nu)(v_j - v_min)(v_max - v_j)`. The parameters of every run are set at the top of the script that
performs it; the defaults of the connectome model are listed in `ConnectomeParameters` (`src/connectome_model.jl`)
and those of the ring in `RingParameters` (`src/ring_model.jl`).

**Interaction function.** In the code the interaction function is written as

```
h(x) = amplitude * ( -sin(x - a) + r sin(2x) ),
```

so that a phase lag alpha in the form `sin(x - alpha)` used in the figure captions corresponds to `a = alpha + pi`.
The values used are

| Used in | h(x) in the figure captions | a | r | amplitude |
|---|---|---|---|---|
| Fig. 1(b),(c); Fig. 2(b) | sin x | pi (3.14159265 in Fig. 2) | 0 | 1 |
| Fig. 2(c) | sin(x - pi/4) | 3.92699082 | 0 | 1 |
| Fig. 2(d) | sin(x - 3pi/8) | 4.31968990 | 0 | 1 |
| Fig. 2(e) | sin(x - pi/2) | 4.71238898 | 0 | 1 |
| Fig. 2(f) | c_1 [sin(x - pi/2) + sin(2x)/3] | 4.71238898 | 0.33333333 | c_1 = 0.86083339 |
| Fig. 1(d); Fig. 2(g); Fig. S1 | c_2 [sin(x - pi/2) + sin 2x] | 4.71238898 | 1 | c_2 = 0.56812611 |
| Fig. 3; Table S1 (ring) | sin(x - pi/2) | 3 pi / 2 | 0 | 1 |

The amplitudes c_1, c_2 normalise h to max h = 1. The numbers are written in the scripts with exactly the digits
used for the published runs.

**Normalisation of the connectome.** The weights are rescaled to an average row sum of 1, so that sigma measures the
total input per node, and the distances to mean 1.

**Initial conditions.** Random initial conditions are theta_i ~ U(0, 2 pi) and v_i ~ N(2.5, 0.1) restricted to
[v_min, v_max], held constant for t <= 0. The time series of Fig. 1, Fig. S1(a) and the long chimera run draw them
after `Random.seed!(1)`; the sweeps use `MersenneTwister(seed)` with seed 1 (Fig. 2) or seeds 1-6 (Fig. S1), so that
seed 1 gives the same initial condition in every grid cell.

## Numerical method

All delay differential equations are integrated with `MethodOfSteps(Tsit5())` from DelayDiffEq.jl, with
abstol = reltol given by

| Computation | Tolerance |
|---|---|
| Fig. 1 time series | 1e-7 |
| Fig. 2 sweeps and neighbourhood checks | 1e-4 |
| Fig. S1 (time series and sweeps), long chimera run | 1e-4 |
| ring (Fig. 3, Table S1) | 1e-8 |

At 1e-4 single nodes of the connectome can make spurious 2 pi phase slips, visible as short dips of R1 in otherwise
synchronous states. The time series of Fig. 1 were therefore computed at 1e-7 (for row (b) a run at 1e-9 gives the
same order parameter). The sweeps
classify the long-time behaviour of thousands of runs and were computed at 1e-4.

## Reproducibility

With the environment of `Manifest.toml` (Julia 1.13.1, Windows 11) the scripts reproduce the included results
exactly. This was checked by rerunning:

* the three Fig. 1 time series and Fig. S1(a): byte-identical CSV files;
* 15 randomly chosen runs of the Fig. 2 sweeps and neighbourhood checks, three of the nine one-parameter sweeps,
  the ring epsilon sweep, and the first 2 x 10^4 time units of the long chimera run: identical values (apart from
  the measured computing times);
* all export scripts: tables identical to those of the published figures, and all four figures pixel-identical
  after compilation;
* `twisted_state_branches.m` (MATLAB R2026a, DDE-Biftool from its GitHub repository): branches and numbers of unstable eigenvalues
  identical to the included results at all 250 points of all six branches.

The green trajectory of Fig. 3 (`figure3/data/frequency_selection_omega_tau_*.csv`) comes from an earlier run with the
same model and parameters; rerunning `simulate_frequency_selection.jl` reproduces it to within 2 x 10^-4, which is
invisible in the figure.

On other platforms or with other package versions, results can differ by rounding errors. Several of the states
studied (chimera, switching) are chaotic, so individual trajectories may then differ after long times, while their
statistics and the classification of the long-time behaviour should not.

## Citation

If you use this code, please cite the Letter above. The connectome data are described in `connectome_data/README.md`.
Simulations use DifferentialEquations.jl (C. Rackauckas and Q. Nie, J. Open Res. Softw. 5, 15 (2017)); stability
computations use DDE-Biftool (J. Sieber, K. Engelborghs, T. Luzyanina, G. Samaey, and D. Roose, arXiv:1406.7144).

## Licence

* Code (`*.jl`, `*.m`, `*.tex`, `*.toml`): [MIT License](LICENSE).
* Simulation results, tables and figures (`data/`, `tables/`, `Figure*.pdf` in the figure folders):
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
* Connectome matrices (`connectome_data/`): redistributed under the
  [WU-Minn HCP Consortium Open Access Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms),
  not under the licences above; see `connectome_data/README.md`.
