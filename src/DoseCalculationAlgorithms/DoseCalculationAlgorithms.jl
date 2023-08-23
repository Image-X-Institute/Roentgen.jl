
"""
    AbstractDoseAlgorithm

Supertype for Dose Algorithms
"""
abstract type AbstractDoseAlgorithm end

"""
    point_dose(p::SVector{3}, beamlet::Beamlet, surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)

Compute the dose at position `p` from `beamlet` using a dose algorithm. 
"""
point_dose

#include("ScaledIsoplaneKernel.jl")
include("FinitePencilBeamKernel.jl")
include("MockKernel.jl")

