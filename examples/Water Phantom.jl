using Pkg; Pkg.activate("."); Pkg.instantiate()
using Roentgen
using Plots, StaticArrays
using SparseArrays, Printf

fieldsize = 100.
SAD = 1000.
SSD = SAD

gantry = GantryPosition(0., 0., SAD);

xb = -0.5*fieldsize:5:0.5*fieldsize
bixels = BixelGrid(xb, xb)
beamlets = Beamlet.(bixels, (gantry,));

calc = FinitePencilBeamKernel("/path/to/kernel-data.jld2")

# Calibrate: 100 MU at maximum with field size of 100 mm and SSD of 1000. mm
MU = 100.
calibrate!(calc, MU, fieldsize, SSD);

surf = PlaneSurface(SSD)

depth = 0.:1.:500
pos = SVector.(0., 0., -depth);

dose = sum(dose_fluence_matrix(Matrix, pos, beamlets, surf, calc); dims=2);

plot(depth, dose, xlabel="depth (mm)", ylabel="Dose (Gy/MU)", legend=false, fmt=:png)

depth = 50.:50:250
x = -300:2:300
pos = SVector.(x, 0., -depth');

dose = sum(dose_fluence_matrix(Matrix, pos, beamlets, surf, calc); dims=2)
dose = reshape(dose, size(pos));

plt = plot(xlabel="x (mm)", ylabel="Dose (Gy/MU)", legendtitle="Depth")
for i in eachindex(depth)
    plot!(plt, x, dose[:, i], label=@sprintf "%0.0fmm" depth[i])
end
plot!(plt; fmt=:png)

Δ = 5.
x = -200.:Δ:200.
y = -200.:Δ:200.
z = -400.:Δ:0.
pos = DoseGrid(x, y, z);

dose = sum(dose_fluence_matrix(Matrix, vec(pos), vec(beamlets), surf, calc); dims=2);
dose = reshape(Array(dose), size(pos));

depth = -z;

k = searchsortedlast(z, -100.)
plt_xy = plot(ylabel="y (mm)", xticks=[])
heatmap!(plt_xy, x, y, dose[:, :, k]', aspect_ratio=1,
         xlim=x[[1,end]], ylim=y[[1, end]])

i = searchsortedlast(x, 0.)
plt_yz = plot(xlabel="z (mm)", yticks=[])
heatmap!(plt_yz, z, y, dose[i, :, :], aspect_ratio=1,
         xlim=z[[1,end]], ylim=y[[1, end]])

j = searchsortedlast(y, 0.)
plt_xz = plot(xlabel="x (mm)", ylabel="depth (mm)")
heatmap!(plt_xz, x, z, dose[:, j, :]', aspect_ratio=1,
         xlim=x[[1,end]], ylim=z[[1, end]])

plot(plt_xy, plt_yz, plt_xz, layout=(2, 2), size=(800, 600); fmt=:png)

# To VTK (Paraview)
write_vtk("water-tank", pos, "dose"=>dose)
# To NRRD (Slicer)
write_nrrd("water-tank.nrrd", pos, dose);
