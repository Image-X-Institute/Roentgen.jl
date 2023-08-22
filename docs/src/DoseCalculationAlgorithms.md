# Dose Calculation Algorithms

Dose calculation algorithms specify how exactly to compute the dose at a given position from a given beamlet.

To create a new dose calculation algorithm, it must be a subtype of [`Roentgen.AbstractDoseAlgorithm`](@ref) and implement the [`Roentgen.point_dose`](@ref) method.
Machine specific parameters should be stored in the struct itself.
This will allow the algorithm to be used with functions such as [`reconstruct_dose`](@ref) and [`dose_fluence_matrix`](@ref)

[`MockKernel`](@ref) is a mock algorithm used as an example and for testing.
Its implementation can be found in [src/DoseCalculationAlgorithms/MockKernel.jl](https://github.com/lmejn/Roentgen.jl/blob/main/src/DoseCalculationAlgorithms/MockKernel.jl), where it implements a simple [`Roentgen.point_dose`](@ref) method.

In addition to [`Roentgen.point_dose`](@ref), it is recommended to implement a [`calibrate!`](@ref) method.
This method calibrates the dose calculation to a given meterset value,
```julia
calibrate!(calc, MU, fieldsize, SAD[, SSD=SAD])
```
Scales the dose such that the maximum dose is 1 Gy for `MU` monitor units, given `fieldsize` and source-axis distance (`SAD`).
The source-surface distance `SSD` can be set if `SSD!=SAD`.

## Finite Pencil Beam Kernel

Implements the [Jelen et al. 2005, "A finite size pencil beam for IMRT dose optimization"](https://dx.doi.org/10.1088/0031-9155/50/8/009) dose calculation kernel.
All references to equations and variables refer to equations in the original paper.

The [`FinitePencilBeamKernel`](@ref) struct stores machine dependent kernel parameters:

- `parameters`: 1D interpolator to obtain $w_1$, $u_x$ and $u_y$ at a given gepth
- `scaling_factor`: 2D interpolator to obtain $A$ at a given depth and tanθ
- `α_depth` and `α_tanθ`: scaling factor to speed up interpolation

The values for $w_1$, $u_x$, $u_y$ and $A$ are typically obtained from Monte-Carlo simulations or water tank data.
Please refer to the original publication for detailed instructions.

[`FinitePencilBeamKernel`](@ref) can be constructed either by directly providing the above parameters,
```julia
calc = FinitePencilBeamKernel(parameters, scaling_factor, α_depth, α_tanθ)
```
or by providing arrays for `parameters`, `scaling_factor`, `depths` and `tanθ`,
```julia
calc = FinitePencilBeamKernel(parameters, scaling_factor, depths, tanθ)
```
These parameters can also be loaded from a [JLD2](https://github.com/JuliaIO/JLD2.jl) filetype (based on the HDF5 standard),
```julia
calc = FinitePencilBeamKernel("path/to/fpbk.jld2")
```
