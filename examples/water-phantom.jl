#=
    Compute the dose on a water phantom for a square field

Compute the dose at various depths in a water phantom, with a 10x10cm² square
field. Uses a plane surface, located at SAD.
=#

using DoseCalculations
using Plots, Meshes, LsqFit

# Set Fieldsize
fieldsize = 100.
ϕg = 0.
θb = 0.
SAD = 1000.

SSD = SAD

# Create Fluence Grid
bixels = [Bixel(0., fieldsize)]

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 280.)

# Calibrate: 100 MU at maximum with field size of 100 mm and SSD of 1000. mm
MU = 100.
calibrate!(calc, MU, fieldsize, SSD)

# Create dose points in IEC BLD coordinates
zp = -SAD-400:1.:-SAD
pos = Vec.(0., 0., zp)

# Create External Surface
surf = PlaneSurface(SAD)

depth = getdepth.(Ref(surf), pos)

# Create dose-fluence matrix
@time dose = dose_fluence_matrix(pos, vec(bixels), surf, calc)

depth = getdepth.(Ref(surf), pos)
begin
    p = plot(xlabel="Depth (mm)", ylabel="Dose (Gy)", ylim=[0, Inf])
    plot!(depth, MU*vec(dose), label="Calculated")
    plot!(depth, MU*norm_depth_dose.(Ref(calc), depth), label="Measured")
    hline!([1.], color=:black, label="")
end
