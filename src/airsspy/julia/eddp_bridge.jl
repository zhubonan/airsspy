#!/usr/bin/env julia

using EDDPotentials
using JSON
using LinearAlgebra
import CellBase

function usage()
    println(stderr,
        "usage: eddp_bridge.jl MODE MODEL INPUT_RES OUTPUT_RES RESULT_JSON " *
        "LABEL METHOD MAX_STEPS F_TOL S_TOL_GPA PRESSURE_GPA RELAX_CELL")
end

if length(ARGS) != 12
    usage()
    exit(2)
end

mode, model_path, input_path, output_path, result_path, label = ARGS[1:6]
method = ARGS[7]
max_steps = parse(Int, ARGS[8])
force_tolerance = parse(Float64, ARGS[9])
stress_tolerance_gpa = parse(Float64, ARGS[10])
pressure_gpa = parse(Float64, ARGS[11])
relax_cell = parse(Bool, ARGS[12])

mode in ("relax", "singlepoint") || error("Unsupported EDDP mode: $mode")

cells = CellBase.read_res_many(input_path)
length(cells) == 1 || error("Expected one structure in $input_path, got $(length(cells))")
cell = deepcopy(only(cells))

calc = load_calculator(model_path)
set_cell!(calc, cell)

converged = true
iterations = 0
fmax = 0.0
smax_gpa = 0.0

if mode == "relax"
    options = RelaxOption(;
        method,
        relax_cell,
        iterations=max_steps,
        convergence="force",
        force_threshold=force_tolerance,
        stress_threshold_gpa=stress_tolerance_gpa,
        external_pressure_gpa=pressure_gpa,
    )
    result = multirelax!(Relax(calc, options); max_iter=max_steps)
    converged = result.converged
    iterations = result.iterations
    fmax = result.fmax
    smax_gpa = result.smax
end

energy = eddp_energy(calc)
forces = eddp_forces(calc)
stress = eddp_stress(calc)
if mode == "singlepoint"
    fmax = maximum(norm.(eachcol(forces)))
    smax_gpa = 160.21766208 * maximum(abs.(stress))
end

write_res(output_path, calc; label)

payload = Dict(
    "mode" => mode,
    "model_path" => abspath(model_path),
    "energy" => energy,
    "forces" => collect.(eachrow(forces)),
    "stress" => collect.(eachrow(stress)),
    "pressure_gpa" => eddp_pressure_gpa(calc),
    "external_pressure_gpa" => pressure_gpa,
    "volume" => volume(get_cell(calc)),
    "converged" => converged,
    "iterations" => iterations,
    "fmax" => fmax,
    "smax_gpa" => smax_gpa,
)

open(result_path, "w") do io
    JSON.print(io, payload)
end
