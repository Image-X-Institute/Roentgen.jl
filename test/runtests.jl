using Test
using DoseCalculations
using HDF5, JLD2
using StaticArrays
using Meshes
using LinearAlgebra, Rotations
using QuadGK, HCubature

@testset "DoseCalculations" begin
    include("utils.jl")
    include("interpolation.jl")
    include("CoordinateSystems.jl")
    include("Gantry.jl")
    include("DosePoints.jl")
    include("DoseFluenceMatrix.jl")
    include("ExternalSurfaces.jl")
    include("meshes.jl")
    include("Bixels.jl")
    include("Beamlet.jl")
    include("Fluence.jl")
    #include("ScaledIsoplaneKernel.jl") # DISABLED, See src/DoseCalculationAlgorithms/ScaledIsoplaneKernel.jl
    include("MultiLeafCollimator.jl")
    include("Jaws.jl")
    include("MultiLeafCollimatorSequence.jl")
    include("TreatmentField.jl")

    # Finite Pencil Beam Kernel
    include("FinitePencilBeamKernel.jl")
end
