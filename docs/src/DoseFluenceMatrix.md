```@setup abc
using Roentgen
```

# Dose-Fluence Matrix

Dose-fluence matrices are precomputed doses for each beamlet used in the treatment.
In dose calculations, it is typical for a broad radiation beam to be split up into a number of beamlets.
These beamlets can then be "turned on/off" depending on whether that beamlet is obscured by a beam-limiting device.

Hence, dose computed through use of a dose-fluence matrix requires four steps:

1. Select the beamlets that are to be used in the treatment
2. Compute the dose at `n` dose positions for `m` beamlets, and store it in a `nxm` matrix, $D$.
3. Compute the fluence and meterset increment for each beamlet, and store in an `m` length vector, $w$.
4. Compute the final dose by multiplying the dose-fluence matrix $D$ by $w$: $\mathrm{dose}=Dw$

An example can be found in the [Dose From Aperture](https://github.com/Image-X-Institute/Roentgen.jl/blob/main/examples/Dose%20from%20Aperture.ipynb) notebook.

## Beamlets

Beamlets have three features:

- An origin, typically where the radiation originates
- A direction, pointing away from the source
- Widths in directions perpendicular to the beamlet direct, usually scaled to the isoplane

Beamlets can be be constructed by combining a `Bixel` and `GantryPosition`,

```@repl abc
bixel = Bixel(0., 0., 1., 1.)
gantry = GantryPosition(0., 0., 1000.)
beamlet = Beamlet(bixel, gantry)
```

This will create a square beamlet of width 1 mm, origin of `[0., 0., 1000.]`, and direction `[0., 0., -1.]`.

Collection of beamlets are created using Julia's broadcast mechanism,

```@repl abc
bixels = BixelGrid(-5.:5.:5., -5.:5.:5.);
beamlets = Beamlet.(bixels, (gantry,))
```

## Creating Dose-Fluence Matrices

Dose-fluence matrices are computed using the [`dose_fluence_matrix`](@ref) and [`dose_fluence_matrix!`](@ref) functions.
These take a dose volume ([Dose Volumes](@ref)), vector of beamlets ([Beamlets](@ref)) and dose calculation algorithm ([Dose Calculation Algorithms](@ref)) and compute the corresponding dose for each position in the volume,
```julia
dose_fluence_matrix(T, vol, beamlets, calc; maxradius=100)
```

By default, dose is not calculated for every element in the dose-fluence matrix.
Dose points that are further away than `maxradius` (scaled to the isoplane) from the beamlet axis are not computed and assumed to be zero.
This speeds up computation (and memory if sparse matrices are used).
`maxradius` is set to 100 mm by default.

[`dose_fluence_matrix`](@ref) creates allocates a new matrix, then computes the matrix values.
The first argument requires the matrix type:

| Type              | Sparse | GPU |
| :---------------- | ------ | --- |
| `Matrix`          | ✖      | ✖   |
| `SparseMatrixCSC` | ✔      | ✖   |
| `CuArray`         | ✖      | ✔   |

Each choice has benefits and downsides.
Using a dense `Matrix` might be the fastest and easiest to manage memory wise, but comes with potentially large memory requirements depending on the number of positions or beamlets.
A sparse `SparseMatrixCSC` only stores non-zero elements of the matrix, and is therefore less memory intensive and may be faster depending on the size of and number of non-zeros in the matrix.

!!! warning
    GPU support is still experimental
If a GPU is available, using `CuArray` will move the computation to the GPU and speed up the computation significantly.
However, all input arguments must be moved to GPU memory, and some types of external surface are not supported on the GPU.
See the [Computing on the GPU](https://github.com/Image-X-Institute/Roentgen.jl/blob/main/examples/Computing%20on%20the%20GPU.ipynb) example for details.

### Performance Intensive Usage

[`dose_fluence_matrix!`](@ref) should be used for more perfomance intensive applications, modifying the supplied matrix.
In this case, the method will specialise depending on the type of the matrix,
```julia
dose_fluence_matrix!(D, vol, beamlets, calc; maxradius=100)
```
where `D` is the matrix.
