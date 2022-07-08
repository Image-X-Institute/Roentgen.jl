#
#   Benchmarks for Calculating Dose-Fluence Matrices
#

using DoseCalculations
using BenchmarkTools, StaticArrays

create_bixels() = bixel_grid([-50., 50.], [-50., 50.], 1.)

function create_mlc_aperture()
    mlc = MultiLeafCollimator(40, 5.)

    mlcx = zeros(2, length(mlc))
    yc = centerposition.(Ref(mlc), 1:length(mlc))
    r = @. √max(0., 30^2 - yc^2)
    @. mlcx[1, :] = -r
    @. mlcx[2, :] = r

    mlcx, mlc
end

function benchmark_dose_fluence_matrix()

    SAD = 1000.

    # Create dose points
    pos = DoseGrid(5., CylinderBounds(100., 100., SVector(0., 0.,  0.)))
    pos_bld = [fixed_to_bld(0., 0., -SAD)(pos[i]) for i in eachindex(pos)]

    # Create an external surface
    surf = PlaneSurface(800.)

    # Choose dose calculation algorithm
    calc = ScaledIsoplaneKernel("src/DoseCalculationAlgorithms/sample-kernel-data/6x/", 25.)

    # Create fluence grid and positions
    bixelgrid = create_bixels()
   
    pos_bld, vec(bixelgrid), surf, calc
end

function benchmark_fluence_from_jaws()
    jaws = Jaws([-32.5, 12.3], [10.2, 40.8])
    bixelgrid = create_bixels()

    bixelgrid, jaws
end

function benchmark_fluence_from_mlc()
    mlcx, mlc = create_mlc_aperture()
    bixelgrid = create_bixels()

    bixelgrid, mlcx, mlc
end

function benchmark_fluence_from_mlc_index()
    mlcx, mlc = create_mlc_aperture()
    bixelgrid = create_bixels()
    index = locate.(Ref(mlc), getindex.(bixelgrid, 2))

    bixelgrid, index, mlcx
end

function benchmark_dose_calculation()
    mlcx, mlc = create_mlc_aperture()

    pos_bld, bixels, surf, calc = benchmark_dose_fluence_matrix()

    D = dose_fluence_matrix(pos_bld, bixels, surf, calc)
    Ψ = fluence(bixels, mlcx, mlc)
    
    D, Ψ
end

function benchmark_dose_reconstruction()
    plan = load_dicom("examples\\sample-data\\RP.zzSPARK_PAT05.PROSTATE.dcm")

    # Resample the plan in 2 MU increments
    plan = resample.(plan, 20.; by=:MU)

    calc = ScaledIsoplaneKernel("src/DoseCalculationAlgorithms/sample-kernel-data/6x/", 25.)

    # Create dose points
    filename = "examples/sample-data/meshes/PTV.stl"
    structure = load_structure_from_ply(filename)
    isocenter = getisocenter(plan[1])

    transform!(structure, patient_to_fixed(isocenter))

    pos = DoseGrid(20., CylinderBounds(100., 100., SVector(0., 0., 0.)))

    # Load surface
    surf = PlaneSurface(800)

    pos, surf, plan, calc
end

# Create benchmark suite
const SUITE = BenchmarkGroup()

SUITE["Dose-Fluence Matrix"] = BenchmarkGroup()
SUITE["Dose-Fluence Matrix"]["creation-parallel"] = @benchmarkable dose_fluence_matrix(args...) setup=( args=benchmark_dose_fluence_matrix() )
SUITE["Dose-Fluence Matrix"]["calculation"] = @benchmarkable args[1]'*args[2] setup=( args=benchmark_dose_calculation())

SUITE["Fluence"] = BenchmarkGroup()
SUITE["Fluence"]["Jaws"] = @benchmarkable fluence(args...) setup=( args=benchmark_fluence_from_jaws() )
SUITE["Fluence"]["MultiLeafCollimator"] = @benchmarkable fluence(args...) setup=( args=benchmark_fluence_from_mlc() )
SUITE["Fluence"]["MultiLeafCollimator Index"] = @benchmarkable fluence(args...) setup=( args=benchmark_fluence_from_mlc_index() )

SUITE["Dose Reconstruction"] = @benchmarkable reconstruct_dose(args...; bixel_size=5., show_progess=false) setup=( args=benchmark_dose_reconstruction() )
