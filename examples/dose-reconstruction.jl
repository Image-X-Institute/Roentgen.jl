#   Calculate the dose of a whole treatment plan
# 
# This script calls quite a high level function, `dose_from_plan`, which
# calculates the dose for all fields in a treatment plan. The dose is returned
# by the function. It also stores the dose at every control point to `output_dir`.

using DoseCalculations, StaticArrays

# Load the plan
plan = load_dicom("path/to/dicom/RP.....dcm")

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 25.)
calibrate!(calc, 100., 100., 1000.)

# Create dose points
grid = DoseGrid(5., CylinderBounds(200., 200., SVector(0., 0., 0.)))

# Load external surface
surf = PlaneSurface(800.)

# Reconstruct Dose
dose = reconstruct_dose(pos, surf, plan, calc)
save("dose", grid, Dict("dose"=>dose))
