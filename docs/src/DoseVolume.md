# Dose Volumes

Dose volumes comprise of two objects, [Dose Positions](@ref) and [External Surfaces](@ref).

Dose positions are essentially lists of points in 3D space where dose is computed.
The external surface is used to compute source-surface distances and depths.

Dose volumes, stored in [`DoseVolume`](@ref), are constructed by first constructing the dose positions and external surfaces.
They can then be used in various computations.

```julia
vol = DoseVolume(pos, surf)
```
