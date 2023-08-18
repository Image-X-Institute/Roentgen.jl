# Dose Volumes

Dose volumes comprise of two objects, dose positions and external surfaces.

Dose positions are essentially lists of points in 3D space where dose is computed.
The external surface is used to compute source-surface distances and depths.

Dose volumes are constructed by first constructing the dose positions and external surfaces.
They can then be used in various computations.

```julia
vol = DoseVolume(pos, surf)
```
