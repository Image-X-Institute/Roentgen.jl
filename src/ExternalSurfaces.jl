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

export ConstantSurface, PlaneSurface, MeshSurface, CylindricalSurface
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

A plane at a constant distance from and normal towards the source

The source-surface distance is stored in `surf.source_surface_distance`
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

Returns `Inf` if no intersection is found. 
"""
getSSD(surf::MeshSurface, pos, src) = getSSD(surf, Point(pos), Point(src))

function getSSD(surf::MeshSurface, pos::Point, src::Point)
    line = Ray(src, pos-src)
    pI, _ = intersect_mesh(line, surf.mesh)
    length(pI)==0 && return Inf
    minimum(norm.(pI .- Ref(src)))
end

write_vtk(filename::String, surf::MeshSurface) = write_vtk(filename, surf.mesh)

#--- CylindricalSurface --------------------------------------------------------

"""
    CylindricalSurface

Surface stored on a cylindrical-polar grid.
"""
struct CylindricalSurface{Ty<:AbstractVector, Tϕ<:AbstractVector, Tdist<:AbstractMatrix, TInterpolation} <: AbstractExternalSurface
    ϕ::Tϕ
    y::Ty
    distance::Tdist
    I::TInterpolation
    function CylindricalSurface(ϕ, y, rho)
        I = LinearInterpolation((ϕ, y), rho)
        new{typeof(ϕ), typeof(y), typeof(rho), typeof(I)}(ϕ, y, rho, I)
    end
end

# Constructors

"""
    CylindricalSurface

Construct from a mesh.
"""
function CylindricalSurface(mesh::SimpleMesh; Δϕ°=2., Δy=2.)
    ϕ = (-180:Δϕ°:180)*π/180

    box = boundingbox(mesh)

    SAD = diagonal(box)

    y₀ = coordinates(minimum(box))[2]
    y₁ = coordinates(maximum(box))[2]
    y = y₀:Δy:y₁
    y = 0.5*(y[2:end]+y[1:end-1])

    rho = zeros(length(ϕ), length(y))

    for j in eachindex(y), i in eachindex(ϕ[1:end-1])
        pos = Point(0., y[j], 0.)
        src = Point(SAD*sin(ϕ[i]), y[j], SAD*cos(ϕ[i]))
        
        line = Ray(src, pos-src)
        pI, _ = intersect_mesh(line, mesh)
        if length(pI)==0
            ρᵢ = Inf
        else
            s = argmin(norm.(pI .- Ref(src)))
            ρᵢ = norm(pI[s]-pos)
        end

        rho[i, j] = ρᵢ
    end
    rho[end, :] .= rho[1, :]

    CylindricalSurface(ϕ, y, rho)
end

# Methods
function distance_to_surface(λ, surf, pos, src)
    r = src + λ*(pos - src)

    ϕ, y = atan(r[1], r[3]), r[2]
    rho = surf.I(ϕ, y)

    idx = SVector(1, 3)
    rho^2 - dot(r[idx], r[idx])
end

function getSSD(surf::CylindricalSurface, pos, src)

    sign(distance_to_surface(0., surf, pos, src)) == sign(distance_to_surface(1., surf, pos, src)) && return Inf

    λ = find_zero(x->distance_to_surface(x, surf, pos, src), (0., 1.), AlefeldPotraShi())
    λ*norm(src-pos)
end

function write_vtk(filename::String, surf::CylindricalSurface)

    x = @. surf.distance*sin(surf.ϕ)
    y = surf.y
    z = @. surf.distance*cos(surf.ϕ)

    xg = reshape(x, size(x)..., 1)
    yg = ones(size(xg)).*y'
    zg = reshape(z, size(z)..., 1)

    vtk = vtk_grid(filename, xg, yg, zg)
    vtk_save(vtk)
    
end

function extent(surf::CylindricalSurface)
    x = @. surf.distance*sin(surf.ϕ)
    y = surf.y
    z = @. surf.distance*cos(surf.ϕ)
    SVector(minimum(x), minimum(y), minimum(z)), SVector(maximum(x), maximum(y), maximum(z))
end

function isinside(surf::CylindricalSurface, pos)
    pos[1]^2+pos[3]^2 == 0. && return true

    ϕ, y = atan(pos[1], pos[3]), pos[2]
    pos[1]^2+pos[3]^2 < surf.I(ϕ, y)^2
end
