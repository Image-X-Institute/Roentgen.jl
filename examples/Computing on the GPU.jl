using Pkg; Pkg.activate("."); Pkg.instantiate()
using Roentgen
using Plots, StaticArrays, CUDA

depth = 0.:2.:400
x = -200.:2.:200.
pos = SVector.(x, 0., -depth');

fieldsize = 100.
SAD = 1000.
SSD = SAD

xb = -0.5*fieldsize:5:0.5*fieldsize
bixels = BixelGrid(xb, xb)

gantry = GantryPosition(0., 0., SAD)
beamlets = Beamlet.(bixels, (gantry,));

calc = FinitePencilBeamKernel("/path/to/kernel-data.jld")
MU = 1. # We'll calibrate such that 1MU=1Gy at max for simplicity
calibrate!(calc, MU, fieldsize, SSD)

surf = PlaneSurface(SSD);

pos = cu(pos);
beamlets = cu(beamlets)
surf = cu(surf)
calc = cu(calc);

CUDA.@sync D = dose_fluence_matrix(CuArray, pos, beamlets, surf, calc);

Ψ = CUDA.fill(1., size(beamlets));
CUDA.@sync dose = D*vec(Ψ);

function plot_dose(dose)
        dose_cpu = Array(dose) # Moves back to the CPU
        heatmap(x, depth, reshape(dose_cpu, size(pos))',
                xlim=x[[1,end]], ylim=depth[[1,end]], aspect_ratio=1,
                xlabel="x (mm)", ylabel="Depth (mm)", clim=(0, 1), fmt=:png)
end
plot_dose(dose)

Ψ = cu(rand(size(beamlets)...))
CUDA.@sync dose = D*vec(Ψ);
plot_dose(dose)
