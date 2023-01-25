#   Calculate the dose of a whole treatment plan
# 
# This script calls quite a high level function, `dose_from_plan`, which
# calculates the dose for all fields in a treatment plan. The dose is returned
# by the function. It also stores the dose at every control point to `output_dir`.

using DoseCalculations

# Load the plan
plan = load_dicom("examples/sample-data/RP.zzSPARK_PAT05.PROSTATE.dcm")

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 25.)
calibrate!(calc, 100., 100., 1000.)

# Load external surface
mesh = load_structure_from_ply("path/to/stl-or-ply")

trans = patient_to_fixed(getisocenter(plan[1]))
surf = CylindricalSurface(transform(mesh, trans))

# Create dose points
pos = DoseGridMasked(5., SurfaceBounds(surf), trans)

# Reconstruct Dose
dose = reconstruct_dose(pos, surf, plan, calc)

# Save to File
write_nrrd("dose.nrrd", pos, dose)
