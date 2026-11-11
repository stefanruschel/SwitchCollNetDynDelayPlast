# setup.jl
#
# Installs the Julia packages of this repository in the exact versions recorded in Manifest.toml
# (DifferentialEquations.jl, DelayDiffEq.jl, MAT.jl and standard libraries). Run once, from this folder:
#
#   julia setup.jl
#
# Afterwards every script is run with  julia --project=.. <script>.jl  from its own folder (see README.md).

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.precompile()
println("Julia environment ready.")
