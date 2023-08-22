# Dose Calculation Algorithms

Dose calculation algorithms specify how exactly to compute the dose at a given position from a given beamlet.

To create a new dose calculation algorithm, it must be a subtype of [`Roentgen.AbstractDoseAlgorithm`](@ref) and implement the [`Roentgen.point_dose`](@ref) method.
Machine specific parameters should be stored in the struct itself.
This will allow the algorithm to be used with functions such as [`reconstruct_dose`](@ref) and [`dose_fluence_matrix`](@ref)


## Mock Kernel

[`MockKernel`](@ref) is a mock algorithm used as an example and for testing.
Its implementation can be found in [src/DoseCalculationAlgorithms/MockKernel.jl](https://github.com/lmejn/Roentgen.jl/blob/main/src/DoseCalculationAlgorithms/MockKernel.jl), where it implements a simple `point_dose` method.
