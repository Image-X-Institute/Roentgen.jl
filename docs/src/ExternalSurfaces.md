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

- Where accuracy is preferred over speed, [Cylindrical Surface](@ref) and [Mesh Surface](@ref) provide accurate source-surface and depth computations.
  Generally Cylindrical Surface is recommended over Mesh Surface.
- In the case of beam commissioning or water-tank doses, [Plane Surface](@ref) is recommended.
- Where speed is preferred, consider using [Linear Surface](@ref)
- For specialised setups, it may be best to implement a custom external surface.

### Plane Surface

A simple plane surface ([`PlaneSurface`](@ref)), positioned at a constant source-surface distance along the central beam axis with normal pointed towards the source.
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

- `rho`: radial distance from the IEC fixed y axis.
- `y`: axial position along the IEC fixed y axis.
- `ϕ`: azimuthal position, equivalent to the gantry angle.

This a faster way of computing mesh intersections, recommended over using a [Mesh Surface](@ref).

[`CylindricalSurface`](@ref) can be constructed by directly supplying the cylindrical-polar surface: `y`, `ϕ` and `rho`,
```julia
ϕ = 0:deg2rad(2.):2π
y = -100.:2.:100.
rho = @. (80+20*sin(ϕ))*(1-(y'/200)^2)
surf = CylindricalSurface(ϕ, y, rho)
```

`CylindricalSurface` can also be constructed using a mesh,
```julia
mesh = load_structure_from_ply("path/to/stl-or-ply")
surf = CylindricalSurface(mesh, Δy)
```
In this example, the bounds in the `y` direction will encompass the whole mesh with spacing of `Δy`.
The axial `y` axis can be provided instead of a spacing, allowing finer control of the extent of the surface.
See [`CylindricalSurface`](@ref) for further details.
!!! note
    `CylindricalSurface` assumes the mesh can be well defined in a cylindrical-polar coordinate system.
    It is up to the user to ensure this, which may require rotating the mesh such.
    See [Coordinate Transformations](@ref) for more information on rotating meshes.

[`CylindricalSurface`](@ref) can be written to the VTK file format for visualisation:
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

[`MeshSurface`](@ref) can be written to the VTK file format for visualisation:
```julia
write_vtk("surface", surf)
```

### Linear Surface

[`LinearSurface`](@ref) is a linear approximation to a general surface.
It computes a set of planes normal to a general surface at the intersection of the line between the iso-center and source positions.

It is recommended to construct from a mesh,
```julia
mesh = load_structure_from_ply("path/to/stl-or-ply")
surf = LinearSurface(mesh)
```
where it computes the plane position and normal for each degree of gantry angle. 
See [`LinearSurface`](@ref) for further details.

### Constant Surface

A constant surface ([`ConstantSurface`](@ref)) returns constant source-surface distance.
Depth is computed as the distance from point to source minus the source-surface distance.
Generally only used for testing purposes.
