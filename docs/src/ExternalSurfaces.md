# External Surfaces

This library provides an interface to provide an external surface, which is used for source-surface distance (SSD) and depth calculations.
Physically, the external surface denotes the boundary between air and the dose absorbing medium (*e.g.* the patient's skin or the surface of a phantom).

All external surfaces are subtypes of `AbstractExternalSurface`, which exposes two methods:

- [`getSSD(surf, p)`](@ref): Returns distance between the surface and the source for a given position in the dose volume (`p`)
- [`getdepth(surf, p)`](@ref): Returns the depth of a given position in the dose volume (`p`)

Both of these methods take an external surface (`surf`) and a dose point position (`pᵢ`). The position is usually a location where dose is computed, and is a three-element vector in the IEC Beam-Limiting-Device (BLD) coordinate system.

All surfaces are assumed to be in the IEC BLD coordinate system. This means that they must undergo coordinate transformation in order to account for gantry and/or collimator rotation.

## Types of External Surface

### Plane Surface

A simple plane surface, positioned at a constant source-surface distance along the central beam axis with normal pointed towards the source. It accounts the fact that SSD increases with off-axis position.

```@repl
using DoseCalculations
surf = PlaneSurface(800.)
getSSD(surf, [0., 0., 810.])
getSSD(surf, [10., 20., 810.])
```

### Mesh Surface

A general surface as defined by a 3D mesh provided by the user.

This surface is created by suppling a 3D `mesh`,
```julia
MeshSurface(mesh)
```
which can be read from a `.ply` file, as described in [Structures](@ref)

This is considered a slow method: for each position in the dose volume, it iterates through every face in the mesh (see [`DoseCalculations.intersect_mesh`](@ref)). A better alternative would be to create an [Isoplane Surface](@ref).


### Isoplane Surface

An isoplane surface collapses a general surface onto a plane centred at the isocentre with normal towards the source. It stores source-surface distances for points on the isoplane. The position of these points are determined by a pre-defined grid.

This makes use of the fact that the distance between source and surface does not change along ray-lines.

The Isoplane surface can be constructed by specifying either the `x` and `y` axes of the grid,
```julia
IsoplaneSurface(x, y, SAD)
```
or by specifying dose volume positions (`pos`), and plane grid spacings `Δx` and `Δy`.
```julia
IsoplaneSurface(pos, Δx, Δy, SAD)
```
In this case, the extent of the grid is computed from dose volume positions.

The SSD values on the plane are computed from a `mesh_surf::DoseCalculations.MeshSurface`,
```julia
compute_SSD!(surf, mesh_surf)
```

### Constant Surface

A constant surface simply returns the provided source-surface distance, regardless of position.

```@repl
using DoseCalculations
surf = ConstantSurface(800.)
getSSD(surf, [0., 0., 810.])
getSSD(surf, [10., 20., 810.])
```

While such a surface is technically possible to implement in reality, this surface is largely used for debug purposes.
