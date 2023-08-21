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
```@setup extsurf
using Roentgen
```
```@repl extsurf
surf = PlaneSurface(800.)
s = [0., 0., 1000.]
p = [10., 20., 0.]
getSSD(surf, p, s)
getdepth(surf, p, s)
```

## Types of External Surface

To enable source-surface and depth computations, a number of concrete external surface types are defined.
These store the surface information and provides the implementation of the `getSSD` and `getdepth` methods.

The choice of which type to use is largely dependent on your use case, whether you prefer  accuracy or speed.

- Where accuracy is prefered over speed, [`CylindricalSurface`](@ref) and [`MeshSurface`](@ref) provide accurate source-surface and depth computations. Generally [`CylindricalSurface`](@ref) is recommended over [`MeshSurface`](@ref).
- In the case of beam commissioning, [`PlaneSurface`](@ref) is recommended.
- For specialised setups, it may be best to implement a custom external surface.

### Plane Surface

A simple plane surface, positioned at a constant source-surface distance along the central beam axis with normal pointed towards the source.
Construct by supplying a source-surface distance, *e.g.* 800mm:
```julia
surf = PlaneSurface(800.)
```

This surface is best used in two scenarios:

1. The external surface is a box and the source does not rotate.
This setup is common in beam commissioning, where dose in computed in a water tank.
2. Accuracy is less important than speed: this method computes source-surface distance or depth much faster than other methods.

### Cylindrical Surface

The [`CylindricalSurface`](@ref) stores the external surface on a [cylindrical-polar](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system) grid.
This coordinate system is aligned to the IEC Fixed coordinate system, defined by the following convention:

- `rho`: Axial distance from the IEC fixed y axis.
- `y`: Signed distance along the IEC fixed y axis.
- `ϕ`: Equivalent to the gantry angle.

This a faster way of computing mesh intersections, recommended over using a [Mesh Surface](@ref).

`CylindricalSurface` can be constructed by directly supplying the cylindrical-polar surface: `y`, `ϕ` and `rho`.
In this example, a surface as defined by `fun` is created:
```julia
ϕ = (-180:2.:180)*π/180
y = -100.:2.:100.
fun(ϕ, y) = (80+20*sin(ϕ))*(1-(y/200)^2)
rho = fun.(ϕ, y')
surf = CylindricalSurface(ϕ, y, rho)
```

`CylindricalSurface` can also be constructed using a mesh:
```julia
mesh = load_structure_from_ply("path/to/stl-or-ply")
surf = CylindricalSurface(mesh; Δϕ°=2., Δy=2.)
```
It will automatically select the bounds in the `y` direction in order to encompass the whole mesh.
Resolutions of the `ϕ` and `y` axes are specified by specifying `Δϕ°` and `Δy`.
!!! note
    `CylindricalSurface` assumes the mesh can be well defined in a cylindrical-polar coordinate system.
    It is up to the user to ensure this, which may require rotating the mesh such.
    See [Coordinate Transformations](@ref) for more information on rotating meshes.

`CylindricalSurface` can be written to the VTK file format for visualisation:
```julia
write_vtk("surface", surf)
```

### Mesh Surface

The [`MeshSurface`](@ref) is defined by a 3D mesh provided by the user.

This surface can be constructed by supplying a mesh:
```julia
mesh = load_structure_from_ply("path/to/stl-or-ply")
surf = MeshSurface(mesh)
```

!!! warning
    This is considered a slow method - the ray-tracing method is inefficient as it iterates through every face in the mesh (see [`Roentgen.intersect_mesh`](@ref)) - and should only be used if the mesh is small, computation time is not important, or a better surface is not available.
    A better alternative would be to create an [Cylindrical Surface](@ref).

`MeshSurface` can be written to the VTK file format for visualisation:
```julia
write_vtk("surface", surf)
```

### Constant Surface

A constant surface returns constant source-surface distance.
Depth is computed as the distance from point to source minus the source-surface distance.
Generally only used for testing purposes.
