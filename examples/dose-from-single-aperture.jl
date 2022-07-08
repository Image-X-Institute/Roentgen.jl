#   Compute the dose on a set of dose points for a single aperture

# This example computes the dose from the first control point in the first
# field of the DICOM plan. It computes on a set of dose points arranged in a
# uniform grid of 2mm spacing in the shape of a cylinder. It uses a simple dose
# calculation kernel, which consists of a gaussian kernel and includes
# beam divergence.

# It plots the fluence map, and save the dose to file in a VTK file format.
# The dose can be visualised using Paraview.

using DoseCalculations
using Plots, Meshes, StaticArrays

# Load the DICOM plan
plan = load_dicom("path/to/dicom/RP.....dcm")
field = plan[1] # Select the first field

ϕg = getϕg(field, 1)
θb = getθb(field, 1)
SAD = getSAD(field)

grid = DoseGrid(5., CylinderBounds(200., 200., SVector(0., 0., 0.)))
pos = [SVector(grid[i]) for i in eachindex(grid)]

# Create Bixels
bixels = bixel_grid(getmlc(field), getjaws(field), 1.)

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 25.)
calibrate!(calc, 100., 100., 1000.)

# Transform positions from IEC Fixed to IEC BLD
trans = fixed_to_bld(ϕg, θb, SAD)
pos_bld = [trans(pos[i]) for i in eachindex(pos)]

# Create External Surface
surf = PlaneSurface(800.)

SSD = getSSD.(Ref(surf), pos_bld)
depth = getdepth.(Ref(surf), pos_bld)

# Create dose-fluence matrix
@time D = dose_fluence_matrix(pos_bld, vec(bixels), surf, calc)

# Compute fluence map from aperture
mlcx = getmlcx(field, 1)

@time Ψ = fluence(bixels, mlcx, getmlc(field)); # Compute the fluence

begin   # Plot fluence map
    p = plot(xlabel="x BLD (mm)", ylabel="y BLD (mm)", aspect_ratio=1)
    plot_bld!(p, bixels, Ψ)
    plot_bld!(p, mlcx, getmlc(field); color=1, fill=false, label="MLC")
    plot_bld!(p, getjaws(field); color=2, label="Jaws")
    axes_lims!(p, getjaws(field))
end

# Compute Dose
@time dose = D'*vec(Ψ)

# Save dose to VTK file format
save("dose", grid, Dict("dose"=>dose,
                        "SSD"=>SSD,
                        "depth"=>depth))

