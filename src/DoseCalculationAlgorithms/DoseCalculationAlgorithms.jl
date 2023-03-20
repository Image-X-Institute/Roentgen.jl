
abstract type AbstractDoseAlgorithm end

#include("ScaledIsoplaneKernel.jl")
include("FinitePencilBeamKernel.jl")

"""
    fluence_grid_length(calc::AbstractDoseAlgorithm)

Return number of bixels in the fluence grid.
"""
fluence_grid_length(calc::AbstractDoseAlgorithm) = length(calc.fluence_grid)

"""
    fluence_grid_size(calc::AbstractDoseAlgorithm)

Return the size of the fluence grid.
"""
fluence_grid_size(calc::AbstractDoseAlgorithm) = size(calc.fluence_grid)
