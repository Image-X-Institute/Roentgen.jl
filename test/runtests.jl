using Test
using DoseCalculations

@testset "DoseCalculations" begin
    include("interpolation.jl")
    include("CoordinateSystems.jl")
    include("ExternalSurfaces.jl")
    include("meshes.jl")
    include("Fluence.jl")
    include("ScaledIsoplaneKernel.jl")
    include("MultiLeafCollimator.jl")
end
