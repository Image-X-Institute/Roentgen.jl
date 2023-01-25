#=
    Compute the dose in a water phantom for a square field

Compute the dose at various depths and off-axis positions in a water phantom
with a 10x10cm² square field. Uses a plane surface, located at SAD.

Computes the dose for three dose profiles:
- The depth-dose curve along the central axis
- The off-axis profile at a depth of 100 mm
- A 2D dose distribution.
=#

using DoseCalculations
using Plots, StaticArrays

#--- Setup Dose Calculation --------------------------------------------------------------------------------------------

# Set Fieldsize
fieldsize = 100.
SAD = 1000.

SSD = SAD

# Create Fluence Grid
bixels = [Bixel(0., fieldsize)]

# Create dose calculation kernel
calc = ScaledIsoplaneKernel("examples/sample-data/dose-kernel/scaled-isoplane-kernel.json", 280.)

# Calibrate: 100 MU at maximum with field size of 100 mm and SSD of 1000. mm
MU = 100.
calibrate!(calc, MU, fieldsize, SSD)

# Create External Surface
surf = PlaneSurface(SAD)

gantry = GantryPosition(0., 0., SAD)
src = DoseCalculations.getposition(gantry)

#--- Depth Dose --------------------------------------------------------------------------------------------------------

# Create dose points
zp = -400:1.:0
pos = SVector.(0., 0., zp)

# Calculate dose
@time depth_dose = dose_fluence_matrix(pos, vec(bixels), gantry, surf, calc)[1, :]

begin
    depth = getdepth.(Ref(surf), pos, Ref(src))
    p = plot(xlabel="Depth (mm)", ylabel="Dose (Gy)", ylim=[0, Inf])
    plot!(depth, MU*vec(depth_dose), label="Calculated")
    plot!(depth, MU*DoseCalculations.norm_depth_dose.(Ref(calc), depth))
    hline!([1.], color=:black, label="")
end

#--- Off-Axis Profile --------------------------------------------------------------------------------------------------

# Create dose points
depth = 100.
xp = -200.:1.:200.
pos = SVector.(xp, 0., -depth)

# Calculate dose
@time offaxis_dose = dose_fluence_matrix(pos, vec(bixels), gantry, surf, calc)[1, :]

# Plot
begin
    p = plot(xlabel="Off-Axis Position (mm)", ylabel="Dose (Gy)", ylim=[0, Inf])
    plot!(xp, MU*offaxis_dose, label="Depth=$depth mm")
end

#--- 2D Dose -----------------------------------------------------------------------------------------------------------

# Create dose points
depth = 0.:5.:400
offaxis_pos = -200.:5.:200.
pos = DoseGrid(offaxis_pos, [0.], -depth)

# Calculate dose
@time dose = MU*dose_fluence_matrix(pos, vec(bixels), gantry, surf, calc)
dose = reshape(Array(dose), size(pos));

# Plot
begin
    p = plot(xlabel="Off-Axis Position (mm)", ylabel="Depth (mm)")
    heatmap!(p, offaxis_pos, depth, dose[:, 1, :]', colorbar_title="Dose (Gy)")
end
