#=
    External surfaces

This file contains stuctures for computing source-surface distances (SSDs) from
external surfaces. External surfaces are structures which encompasses all dose
points, and are hence used to compute the depth of the dose point within the
surface.

It current contains two types of external surfaces:

 - ConstantSurface: Returns a constant SSD regardless of position of the dose
                    point or orientation of the treatment beam.

 - PlaneSurface: Returns the distance from the source to a plane, but does not
                 account for the orientation of the treatment beam.
=#

export ConstantSurface, PlaneSurface, MeshSurface, VariablePlaneSurface
export getdepth, getSSD, compute_SSD!

#--- AbstractExternalSurface --------------------------------------------------
"""
    AbstractExternalSurface

Supertype for external surfaces.
"""
abstract type AbstractExternalSurface end

"""
    getSSD(surf::AbstractExternalSurface, pos, src)

Get the Source-Surface Distance (SSD) for position `pos` to the radiation source `src`.
"""
getSSD(surf::AbstractExternalSurface, pos, src)

"""
    getdepth(surf::AbstractExternalSurface, pos, src)

Get the depth of the position `pos` below the surface from the radiation source `src`.

Computes the depth by subtracting 
"""
getdepth(surf::AbstractExternalSurface, pos, src) = norm(pos - src) - getSSD(surf, pos, src)

#--- ConstantSurface ----------------------------------------------------------

"""
    ConstantSurface

An external surface with a constant source surface distance, regardless of the
position. This surface is largely used for testing purposes as its inherently
unphysical.
"""
struct ConstantSurface{T} <: AbstractExternalSurface
    source_surface_distance::T
end

"""
    getSSD(calc::ConstantSurface, pos, src)

When applied to a `ConstantSurface`, it returns a constant source surface distance,
regardless of `pos`.
"""
getSSD(surf::ConstantSurface, pos, src) = surf.source_surface_distance


#--- PlaneSurface -------------------------------------------------------------

"""
    PlaneSurface

A planar external surface at a constant distance from the source.

It assumes the external surface is a plane located at a distance of `surf.source_surface_distance`
away from the source with normal from isocenter to source.
"""
struct PlaneSurface{T} <: AbstractExternalSurface
    source_surface_distance::T
end

"""
    hypotenuse(a, b)

Compute the hypotenuse of triangle
"""
hypotenuse(a, b) = √(dot(a,a)*dot(b,b))/dot(a,b)

"""
    getSSD(surf::PlaneSurface, pos, src)

When applied to a `PlaneSurface`, it returns the distance to the plane.
"""
getSSD(surf::PlaneSurface, pos, src) = surf.source_surface_distance*hypotenuse(src, src - pos)

getSSD(surf::PlaneSurface, pos::Point, src) = getSSD(surf, coordinates(pos), src)

#--- MeshSurface --------------------------------------------------------------

"""
    MeshSurface

An external surface defined by a 3D mesh.
"""
struct MeshSurface{T} <: AbstractExternalSurface
    mesh::Mesh{3, T}
end

"""
    getSSD(surf::MeshSurface, pos, src)

When applied to a `MeshSurface`, it returns the smallest distance to the mesh.
"""
getSSD(surf::MeshSurface, pos, src) = getSSD(surf, Point(pos), Point(src))

function getSSD(surf::MeshSurface, pos::Point, src::Point)
    line = Ray(src, pos-src)
    pI, _ = intersect_mesh(line, surf.mesh)
    length(pI)==0 && return 0.
    minimum(norm.(pI .- Ref(src)))
end

#--- VariablePlaneSurface ------------------------------------------------------

"""
    VariablePlaneSurface

A planar external surface at a variable distance from the isocenter.
"""
struct VariablePlaneSurface{Tinterp} <: AbstractExternalSurface
    distance::Tinterp
end

# Constructors
VariablePlaneSurface(ϕ, distance) = VariablePlaneSurface(LinearInterpolation(ϕ, distance))

# Methods
interpolate(surf::VariablePlaneSurface, src) = surf.distance(atan(src[3], src[1]))
getSSD(surf::VariablePlaneSurface, pos, src) = interpolate(surf, src)*hypotenuse(src, src - pos)
