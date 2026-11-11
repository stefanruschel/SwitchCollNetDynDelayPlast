# export_figure2.jl
#
# Turns the sweeps of run_figure2_sweeps.jl into the plain-text tables that Figure2.tex (pgfplots) draws.
# No simulation is done here.
#
# Usage (from this folder):
#   julia --project=.. export_figure2.jl        then        lualatex Figure2.tex
#
# Inputs:  data/sweep_figure2b.csv ... data/sweep_figure2g.csv
# Outputs: tables/map_<panel>.dat      cell centres (x = log10 sigma, y = nu) and cell colours
#          tables/sw_hatch_<panel>.dat, sw_bnd_<panel>.dat   hatching and boundary of the switching region
#          tables/ch_hatch_<panel>.dat, ch_bnd_<panel>.dat   hatching and boundary of the chimera-like region
#          ("nan nan" lines separate line segments)
#          tables/key.dat (colour key), tables/params.tex (axis limits)
#
# Encoding of a grid cell (the caption of Fig. 2):
#   colour        bilinear interpolation over (<R1>, <R2>) between the corner colours
#                 (0, 0) 50 % grey (incoherent), (0, 1) blue (two-cluster state), (1, 0) and (1, 1) white (synchronous)
#   switching     at least two collapses of R1 or of R2 (rise above 0.9, then fall below 0.1); cyan, hatch '/'
#   chimera-like  0.1 < fraction of frequency-locked nodes < 0.9 and 0.2 < <R1> < 0.9; magenta, hatch '\'
# With several seeds per cell (not used in Fig. 2) a cell is marked if at least half of its runs qualify.

using DelimitedFiles, Statistics, Printf

const TABLES = joinpath(@__DIR__, "tables")
const PANELS = ["b", "c", "d", "e", "f", "g"]
const HATCH_LINES = 10          # hatch density (lines per panel width)

# corner colours of the key (RGB)
const COLOUR_00 = (0.5, 0.5, 0.5)     # <R1> = 0, <R2> = 0
const COLOUR_01 = (0.0, 0.0, 1.0)     # <R1> = 0, <R2> = 1
const COLOUR_11 = (1.0, 1.0, 1.0)     # <R1> = 1, <R2> = 1
const COLOUR_10 = (1.0, 1.0, 1.0)     # <R1> = 1, <R2> = 0

"Colour of a cell with time-averaged order parameters R1, R2 (bilinear interpolation between the corners)."
function cell_colour(R1, R2)
    x = clamp(R1, 0, 1)
    y = clamp(R2, 0, 1)
    weights = ((1 - x) * (1 - y), (1 - x) * y, x * y, x * (1 - y))
    corners = (COLOUR_00, COLOUR_01, COLOUR_11, COLOUR_10)
    return ntuple(k -> weights[1] * corners[1][k] + weights[2] * corners[2][k] + weights[3] * corners[3][k] + weights[4] * corners[4][k], 3)
end

rgb_string(colour) = @sprintf("rgb=%.3f,%.3f,%.3f", colour...)

"Switching: at least two strict collapses of R1 or of R2."
is_switching(collapses_R1, collapses_R2) = collapses_R1 >= 2 || collapses_R2 >= 2

"Chimera-like: a frequency-locked and a drifting group coexist and the order parameter is intermediate."
is_chimera_like(fraction_locked, R1_mean) = 0.1 < fraction_locked < 0.9 && 0.2 < R1_mean < 0.9

"""
    hatching_and_boundary(x, y, dx, dy, region; direction = 1)

Line segments that hatch the cells of the Boolean matrix `region` (rows = y, columns = x) at 45 degrees in panel
coordinates (direction = +1 for '/', -1 for '\'), and the segments of the boundary of the region.
"""
function hatching_and_boundary(x, y, dx, dy, region; direction = 1)
    nx, ny = length(x), length(y)
    width = x[end] - x[1] + 2dx
    height = y[end] - y[1] + 2dy
    x0, y0 = x[1] - dx, y[1] - dy
    hatching = Float64[]
    boundary = Float64[]
    for offset in range(-height, 2height; length = 3 * HATCH_LINES + 1), i in 1:ny, j in 1:nx
        region[i, j] || continue
        left, right, bottom, top = x[j] - dx, x[j] + dx, y[i] - dy, y[i] + dy
        # the hatch line is y = y0 + offset + direction * s * height, x = x0 + s * width; clip s to the cell
        s_lo = (left - x0) / width
        s_hi = (right - x0) / width
        if direction == 1
            s_lo = max(s_lo, (bottom - y0 - offset) / height)
            s_hi = min(s_hi, (top - y0 - offset) / height)
        else
            s_lo = max(s_lo, (y0 + offset - top) / height)
            s_hi = min(s_hi, (y0 + offset - bottom) / height)
        end
        s_lo < s_hi || continue
        append!(hatching, (x0 + s_lo * width, y0 + offset + direction * s_lo * height,
                           x0 + s_hi * width, y0 + offset + direction * s_hi * height, NaN, NaN))
    end
    inside(i, j) = (1 <= i <= ny && 1 <= j <= nx) ? region[i, j] : false
    for i in 1:ny, j in 1:nx
        region[i, j] || continue
        inside(i, j - 1) || append!(boundary, (x[j] - dx, y[i] - dy, x[j] - dx, y[i] + dy, NaN, NaN))
        inside(i, j + 1) || append!(boundary, (x[j] + dx, y[i] - dy, x[j] + dx, y[i] + dy, NaN, NaN))
        inside(i - 1, j) || append!(boundary, (x[j] - dx, y[i] - dy, x[j] + dx, y[i] - dy, NaN, NaN))
        inside(i + 1, j) || append!(boundary, (x[j] - dx, y[i] + dy, x[j] + dx, y[i] + dy, NaN, NaN))
    end
    return permutedims(reshape(hatching, 2, :)), permutedims(reshape(boundary, 2, :))
end

"Writes x y pairs; a NaN row becomes 'nan nan' (a break between line segments)."
function write_segments(name, points)
    open(joinpath(TABLES, name), "w") do io
        println(io, "x y")
        if size(points, 1) == 0
            println(io, "nan nan")
            return
        end
        for point in eachrow(points)
            if isnan(point[1])
                println(io, "nan nan")
            else
                @printf(io, "%.5f %.5f\n", point[1], point[2])
            end
        end
    end
end

function main()
    mkpath(TABLES)
    limits = nothing
    for panel in PANELS
        file = joinpath(@__DIR__, "data", "sweep_figure2$(panel).csv")
        table, header = readdlm(file, ','; header = true)
        header = vec(header)
        column(name) = table[:, findfirst(==(name), header)]

        sigmas = sort(unique(column("sigma")))
        nus = sort(unique(column("nu")))
        x = log10.(sigmas)
        y = collect(nus)
        nx, ny = length(x), length(y)
        dx = (x[2] - x[1]) / 2
        dy = (y[2] - y[1]) / 2
        switching = falses(ny, nx)
        chimera_like = falses(ny, nx)

        open(joinpath(TABLES, "map_$(panel).dat"), "w") do io
            println(io, "x y c")
            for i in 1:ny, j in 1:nx                       # x varies fastest (pgfplots: mesh/ordering = x varies)
                runs = (column("sigma") .≈ sigmas[j]) .& (column("nu") .≈ nus[i])
                R1 = mean(column("R1mean")[runs])
                R2 = mean(column("R2mean")[runs])
                switching[i, j] = mean(is_switching.(column("events_R1_strict")[runs], column("events_R2_strict")[runs])) >= 0.5
                chimera_like[i, j] = mean(is_chimera_like.(column("frac_locked")[runs], column("R1mean")[runs])) >= 0.5
                @printf(io, "%.5f %.5f %s\n", x[j], y[i], rgb_string(cell_colour(R1, R2)))
            end
        end
        for (name, region, direction) in (("sw", switching, 1), ("ch", chimera_like, -1))
            hatching, boundary = hatching_and_boundary(x, y, dx, dy, region; direction = direction)
            write_segments("$(name)_hatch_$(panel).dat", hatching)
            write_segments("$(name)_bnd_$(panel).dat", boundary)
        end
        println("($panel) switching: $(count(switching)) cells, chimera-like: $(count(chimera_like)) cells, of $(nx * ny)")
        limits = (x[1] - dx, x[end] + dx, y[1] - dy, y[end] + dy, nx, ny)
    end

    open(joinpath(TABLES, "key.dat"), "w") do io
        println(io, "x y c")
        n = 48
        for i in 1:n, j in 1:n
            R1 = (j - 0.5) / n
            R2 = (i - 0.5) / n
            @printf(io, "%.5f %.5f %s\n", R1, R2, rgb_string(cell_colour(R1, R2)))
        end
    end
    open(joinpath(TABLES, "params.tex"), "w") do io
        @printf(io, "\\def\\Xmin{%.5f}\\def\\Xmax{%.5f}\\def\\Ymin{%.5f}\\def\\Ymax{%.5f}\\def\\NX{%d}\\def\\NY{%d}\n", limits...)
        @printf(io, "\\def\\XtA{%.5f}\\def\\XtB{%.5f}\\def\\XtC{%.5f}\n", log10(0.3), log10(1.0), log10(3.0))
    end
    println("wrote tables/*")
end

main()
