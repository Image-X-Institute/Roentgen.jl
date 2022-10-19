#   Compute the dose on a set of dose points for a single aperture

# This example computes the dose from the first control point in the first
# field of the DICOM plan. It computes on a set of dose points arranged in a
# uniform grid of 2mm spacing in the shape of a cylinder. It uses a simple dose
# calculation kernel, which consists of a gaussian kernel and includes
# beam divergence.

# It plots the fluence map, and save the dose to file in a VTK file format.
# The dose can be visualised using Paraview.

using DoseCalculations, StaticArrays

# Load the DICOM plan
plan = load_dicom("path/to/dicom/RP.....dcm")
field = plan[1] # Select the first field
beams = resample(field, 2.) # Resample the first field for better accuracy
controlpoint = beams[1] # Select the first control point

pos = DoseGrid(5., CylinderBounds(200., 200., SVector(0., 0., 0.)))

# Create Bixels
bixels = bixel_grid(getmlc(controlpoint), getjaws(controlpoint), 1.)

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 25.)
calibrate!(calc, 100., 100., 1000.)

# Create External Surface
surf = PlaneSurface(800.)

# Create dose-fluence matrix
gantry = getgantry(controlpoint)
@time D = dose_fluence_matrix(pos, vec(bixels), gantry, surf, calc)

# Compute fluence map from aperture
Ψ = fluence(bixels, getmlc(controlpoint)); # Compute the fluence

# Plot fluence map
begin
    p = plot_bld(bixels, Ψ, xlabel="x BLD (mm)", ylabel="y BLD (mm)", aspect_ratio=1)
    plot_bld!(p, getmlc(controlpoint); fill=false, color=1, label="MLC")
    plot_bld!(p, getjaws(controlpoint); color=2, label="Jaws")
    axes_lims!(p, getjaws(controlpoint), pad=10.)
end

# Compute Dose
ΔMU = getΔMU(controlpoint)
dose = ΔMU*D'*vec(Ψ)

# Save dose to VTK file format
save("dose", pos, Dict("dose"=>dose))
