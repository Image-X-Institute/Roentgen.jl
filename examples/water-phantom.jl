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
using SparseArrays, Printf

#--- Setup Dose Calculation --------------------------------------------------------------------------------------------

# Set Fieldsize
fieldsize = 100.
SAD = 1000.

SSD = SAD

# Create Fluence Grid
xb = -0.5*fieldsize:5:0.5*fieldsize
bixels = bixel_grid(xb, xb)

# Create dose calculation kernel
calc = FinitePencilBeamKernel("path/to/kernel/file.jld")

# Calibrate: 100 MU at maximum with field size of 100 mm and SSD of 1000. mm
MU = 100.
calibrate!(calc, MU, fieldsize, SSD)

# Create External Surface
surf = PlaneSurface(SSD)

gantry = GantryPosition(0., 0., SAD)
beamlets = Beamlet.(bixels, (gantry,))

#--- Depth Dose --------------------------------------------------------------------------------------------------------

# Create dose points
depth = 0.:1.:500
pos = SVector.(0., 0., -depth)

# Calculate dose
dose = sum(dose_fluence_matrix(Matrix, pos, beamlets, surf, calc); dims=2)

plot(depth, dose, xlabel="Depth (mm)", ylabel="Dose (Gy/MU)", legend=false)

#--- Off-Axis Profile --------------------------------------------------------------------------------------------------

# Create dose points
depth = 50.:50:250
x = -300:2:300
pos = SVector.(x, 0., -depth')

# Calculate dose
dose = sum(dose_fluence_matrix(Matrix, pos, beamlets, surf, calc); dims=2)
dose = reshape(dose, size(pos))

# Plot
begin

    plt = plot(xlabel="x (mm)", ylabel="Dose (Gy/MU)", legendtitle="Depth")
    for i in eachindex(depth)
        plot!(plt, x, dose[:, i], label=@sprintf "%0.0fmm" depth[i])
    end
    plt
end

#--- 2D Dose -----------------------------------------------------------------------------------------------------------

# Create dose points
depth = 0.:2.:400
x = -200.:2.:200.
pos = SVector.(x, 0., -depth')

# Calculate dose
dose = sum(dose_fluence_matrix(Matrix, pos, beamlets, surf, calc); dims=2)
dose = reshape(dose, size(pos))

# Plot
begin
    plt = plot(xlabel="x (mm)", ylabel="Depth (mm)", aspect_ratio=1,
               xlim=x[[1, end]], ylim=depth[[1, end]])
    heatmap!(plt, x, depth, dose', colorbar_title="Dose (Gy/MU)")
end

#--- 3D Dose -----------------------------------------------------------------------------------------------------------

# Create dose points
Δ = 2.
depth = 0.:Δ:400
x = -200.:Δ:200.
y = -200.:Δ:200.
z = -depth
pos = DoseGrid(x, y, z)

# Calculate dose
dose = sum(dose_fluence_matrix(Matrix, pos, vec(beamlets), surf, calc); dims=2);
dose = reshape(Array(dose), size(pos));

write_vtk("water", pos, "dose"=>dose)

# Plot Slices
begin
    k = searchsortedlast(depth, 100.)
    plt_xy = plot(ylabel="y (mm)", xticks=[])
    heatmap!(plt_xy, x, y, dose[:, :, k]', aspect_ratio=1,
             xlim=x[[1,end]], ylim=y[[1, end]])
    
    i = searchsortedlast(x, 0.)
    plt_yz = plot(xlabel="depth (mm)", yticks=[])
    heatmap!(plt_yz, depth, y, dose[i, :, :], aspect_ratio=1,
             xlim=depth[[1,end]], ylim=y[[1, end]])

    j = searchsortedlast(y, 0.)
    plt_xz = plot(xlabel="x (mm)", ylabel="depth (mm)")
    heatmap!(plt_xz, x, depth, dose[:, j, :]', aspect_ratio=1,
             xlim=x[[1,end]], ylim=depth[[1, end]])

    plot(plt_xy, plt_yz, plt_xz; layout=(2, 2), size=(800, 600))
end
