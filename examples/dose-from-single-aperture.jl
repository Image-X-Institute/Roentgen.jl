#   Compute the dose on a set of dose points for a single aperture

# This example computes the dose from the first control point in the first
# field of the DICOM plan. It computes on a set of dose points arranged in a
# uniform grid of 2mm spacing in the shape of a cylinder. It uses a simple dose
# calculation kernel, which consists of a gaussian kernel and includes
# beam divergence.

# It plots the fluence map, and save the dose to file in a VTK file format.
# The dose can be visualised using Paraview.

using DoseCalculations

# Load the DICOM plan
plan = load_dicom("path/to/dicom/RP.....dcm")
field = plan[1] # Select the first field
controlpoint = field[1] # Select the first control point

# External Surface
mesh = load_structure_from_ply("path/to/stl-or-ply")
trans = patient_to_fixed(getisocenter(controlpoint))

surf = CylindricalSurface(transform(mesh, trans))

# Dose Positions
pos = DoseGridMasked(5., SurfaceBounds(surf), trans)
pos_fixed = trans.(pos)

# Create Bixels
bixels = bixel_grid(getmlc(controlpoint), getjaws(controlpoint), 1.)

# Compute fluence map from aperture
Ψ = fluence(bixels, getmlc(controlpoint)); # Compute the fluence
ΔMU = getΔMU(controlpoint)

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 25.)
calibrate!(calc, 100., 100., 1000.)

# Create dose-fluence matrix
gantry = getgantry(controlpoint)
@time D = dose_fluence_matrix(pos_fixed, vec(bixels), gantry, surf, calc)
dose = ΔMU*D'*vec(Ψ)

# Save dose to VTK file format
write_vtk("dose", pos, "dose"=>dose)
