using Test
using DoseCalculations
using HDF5
using StaticArrays
using Meshes
using LinearAlgebra, Rotations
using QuadGK, HCubature

@testset "DoseCalculations" begin
    include("utils.jl")
    include("interpolation.jl")
    include("CoordinateSystems.jl")
    include("DosePoints.jl")
    include("ExternalSurfaces.jl")
    include("meshes.jl")
    include("Fluence.jl")
    include("ScaledIsoplaneKernel.jl")
    include("MultiLeafCollimator.jl")
    include("MultiLeafCollimatorSequence.jl")
end
