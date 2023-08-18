# using Pkg; Pkg.activate(".")
using Roentgen
using CoordinateTransformations, Rotations

plan = load_dicom("examples/sample-data/RP.zzSPARK_PAT05.PROSTATE.dcm")

mesh = load_structure_from_ply("examples/sample-data/meshes/body.stl")
trans = patient_to_fixed(getisocenter(plan[1]))


calc = FinitePencilBeamKernel("examples/sample-data/fpbk-fieldsize_120mm.jld2")
calibrate!(calc, 100., 100., 1000.);

# surf = PlaneSurface(800.)
cylsurf = CylindricalSurface(transform(mesh, trans), 5.)

surf = LinearSurface(transform(mesh, trans))
Δ = 2.5
pos = DoseGridMasked(Δ, SurfaceBounds(cylsurf), trans);
pos_fixed = trans.(pos);

vol = Roentgen.DoseVolume(pos_fixed, cylsurf)


# pos = DoseGridMasked(2.5, CylinderBounds(200., 200.), LinearMap(RotX(π/2)));

dose_linsurf = reconstruct_dose(pos, surf, plan, calc)
dose_cylsurf = reconstruct_dose(pos, cylsurf, plan, calc)

write_vtk("dose", pos, "dose_linsurf"=>dose_linsurf)#, "dose_cylsurf"=>dose_cylsurf)

write_nrrd("dose_linsurf.nrrd", pos, dose_linsurf)
write_nrrd("dose_cylsurf", pos, dose_cylsurf)

# Generate beamlets

begin
    beamlets = Beamlet{Float64}[]
    ΔMU = Float64[]
    for field in plan, beam in field
        mlc = getmlc(beam)
        jaws = getjaws(beam)
        gantry = getgantry(beam)
        bixels = bixels_from_bld(mlc, jaws; Δx=5.)
        append!(beamlets, Beamlet.(bixels, Ref(gantry)))
        append!(ΔMU, fill(getΔMU(beam), length(bixels)))
    end
    length(beamlets), length(ΔMU)
end

using CUDA
CUDA.allowscalar(false)

cu_vol = cu(vol);
cu_calc = cu(calc);

cu_beamlets = cu(beamlets)

cu_ΔMU = cu(ΔMU);


function gpu_point_dose!(dose, ΔMU, vol, beamlets, calc, maxradius)

    pos = Roentgen.getpositions(vol)
    surf = Roentgen.getsurface(vol)

    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    sx = gridDim().x * blockDim().x

    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    sy = gridDim().y * blockDim().y

    for i = ix:sx:length(pos)
        val = 0.f0
        for j=iy:sy:length(beamlets)
            @inbounds val += ΔMU[j]*Roentgen.point_dose(pos[i],
                beamlets[j], surf, calc, maxradius)
        end
        CUDA.@atomic dose[i] += val
    end

    return nothing
end

N, M = length(pos), length(beamlets)
cu_dose = CUDA.fill(0f0, N);

N*M

kernel = @cuda launch=false gpu_point_dose!(cu_dose, cu_ΔMU, cu_vol, cu_beamlets, cu_calc, 100.f0)
config = launch_configuration(kernel.fun)
threads = min(N, config.threads)
blocks = cld(N, threads)

threads*blocks, N
threads*blocks>=N

@time CUDA.@sync begin
    CUDA.fill!(cu_dose, 0f0)
    kernel(cu_dose, cu_ΔMU, cu_vol, cu_beamlets, cu_calc, 100.f0; threads, blocks)
end

threads = min.((N, M), threads)
blocks = cld.((N, M), threads)

# Each block has
prod(threads)*prod(blocks), N*M
prod(threads)*prod(blocks)>=N*M

@time CUDA.@sync begin
    CUDA.fill!(cu_dose, 0f0)
    kernel(cu_dose, cu_ΔMU, cu_vol, cu_beamlets, cu_calc, 100.f0; threads, blocks)
end

write_vtk("dose", pos, "dose"=>vec(Array(cu_dose)))


84*384
threads*blocks



@time CUDA.@sync begin
    CUDA.fill!(cu_dose, 0f0)
    @cuda kernel(cu_dose, cu_ΔMU, cu_pos, cu_beamlets, cu_surf, cu_calc, 100.f0; threads, blocks)
end

# dose = reconstruct_dose(pos, surf, plan, calc; show_progress=false);

write_vtk("dose", pos, "dose"=>vec(Array(cu_dose)))
