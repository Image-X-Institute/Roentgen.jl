# External Surfaces

This library provides an interface to provide an external surface, which is used for source-surface distance (SSD) and depth calculations.
Physically, the external surface denotes the boundary between air and the dose absorbing medium (*e.g.* the patient's skin or the surface of a phantom).

All external surfaces exposes two methods,

- [`getSSD(surf, p, s)`](@ref): Returns the source-surface distance (SSD)
- [`getdepth(surf, p, s)`](@ref): Returns the depth below the source-surface intersection

![external_surface](assets/external-surface.svg)

Both of these methods take an external surface (`surf`) and two position vectors, `p` and `s`.
The position `p` is usually a location where dose is computed, and is a three-element vector.
The position `s` is the position of the radiation source.
Both positions are in the IEC Fixed coordinate system.

For example:
```@repl
using DoseCalculations
surf = PlaneSurface(800.)
s = [0., 0., 1000.]
p = [10., 20., 0.]
getSSD(surf, p, s)
getdepth(surf, p, s)
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
