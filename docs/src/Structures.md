# Structures

Dose volumes often contain structures of interest (*e.g.* the PTV or an OAR) where one may want to compute dose, or determine various metrics.

This section describes how to use structures to:
* Load structures from file.
* Create dose points within a given structure 
* Tag dose points based on which structure they are in.

## Creating Structures

One can load a mesh file using [`load_structure_from_ply`](@ref). Currently, only `.ply` files are supported.

```julia
structure = load_structure_from_ply("mesh.ply")
```
This returns a `Meshes.SimpleMesh`, as provided by [Meshes.jl](https://juliageometry.github.io/Meshes.jl/stable/meshes.html).

Currently, loading structures direct from a DICOM RT Structure Set is not supported. Structures stored in DICOM files are not in a "friendly" mesh format, and usually require some modification before they can be used. There are many tools available which convert DICOM structures into a suitable mesh file (*e.g.* [3D Slicer](https://www.slicer.org/)).

## Coordinate Transformations

One can apply general transformations to any structure by calling the [`transform`](@ref) or [`transform!`](@ref) function:
```julia
T = fixed_to_bld(0., 0., 1000.)
transform(structure, T)
```
`T` is a general transformation as provided by [CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl). In this example, the helper function [`fixed_to_bld`](@ref) is used to create a transformation from the IEC Fixed to IEC BLD coordinate system. See [Coordinate Transformations](@ref) for more details.
