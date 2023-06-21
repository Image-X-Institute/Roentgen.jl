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
export LinearSurface
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
getSSD(surf::AbstractExternalSurface, pos, src) = norm(pos-src) - getdepth(surf, pos, src)

"""
    getdepth(surf::AbstractExternalSurface, pos, src)

Get the depth of the position `pos` below the surface from the radiation source `src`.

Computes the depth by subtracting 
"""
getdepth(surf::AbstractExternalSurface, pos, src)

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
struct MeshSurface{TMesh, TBox} <: AbstractExternalSurface
    mesh::TMesh
    boxes::Vector{TBox}

    function MeshSurface(mesh::Partition)
        boxes = boundingbox.(mesh)
        new{typeof(mesh), eltype(boxes)}(mesh, boxes)
    end
end

function _boxwidth(mesh)
    box = boundingbox(mesh)
    pmin, pmax = extrema(box)
    pmax - pmin
end

function MeshSurface(mesh::SimpleMesh, boxwidths::AbstractVector{T}) where T<:Real
    part = BlockPartition(boxwidths)
    blockmesh = partition(mesh, part)
    MeshSurface(blockmesh)
end

function MeshSurface(mesh::SimpleMesh, n::Union{Int, AbstractVector{Int}}=8)
    widths = SVector(_boxwidth(mesh))
    MeshSurface(mesh, widths./n)
end

function intersection_points(surf::MeshSurface, pos::T, src::T) where T<:Point
    intersect_mesh(Segment(pos, src), surf.mesh, surf.boxes)
end

function getSSD(surf::MeshSurface, pos, src)
    pt = closest_intersection(pos, src, surf.mesh, surf.boxes)
    pt===nothing && return Inf
    norm(pt-src)
end

function getdepth(surf::MeshSurface, pos::T, src::T) where T<:Point
    pts = intersection_points(surf, pos, src)
    length(pts)==0 && return NaN

    push!(pts, pos)
    sum(@. norm(pts[1:2:end]-pts[2:2:end]))
end
getdepth(surf::MeshSurface, pos, src) = getdepth(surf, Point(pos), Point(src))

#--- Linear Surface ------------------------------------------------------------

struct Plane{T}
    p::SVector{3, T}
    n::SVector{3, T}
end

function intersection_point(plane::Plane, p1::SVector{3, T}, p2::SVector{3, T}) where T<:Real
    v = p2-p1
    nv = dot(plane.n, v)
    isapprox(nv, zero(T), atol=T(1e-16)) && return nothing
    λ = dot(plane.n, plane.p - p1)/nv
    p1 + λ*v
end

struct LinearSurface{T<:AbstractInterpolation}
    params::T

    function LinearSurface(params)
        @assert length(params)==361 "Distance must be supplied at every degree"
        I = interpolate(params, BSpline(Linear()))
        new{typeof(I)}(I)
    end
end

LinearSurface(p, n) = LinearSurface(vcat.(p, n))

function LinearSurface(ϕg, p, n)
    params = vcat.(p, n)
    I = linear_interpolation(ϕg, params)
    dist = I.(deg2rad.(0:360))
    LinearSurface(dist)
end

"""
    LinearSurface(mesh[, SAD=1000, ΔΘ=deg2rad(1)])

Construct a LinearSurface from a mesh.

Computes a set of planes parallel to the surface of the mesh.
"""
function LinearSurface(mesh::SimpleMesh{3, T}; SAD=T(1000.), ΔΘ=deg2rad(1)) where {T<:Real}
    N = 361
    ϕg = 2π*range(0, 1, length=N)
    n = Vector{SVector{3, T}}(undef, N)
    p = Vector{SVector{3, T}}(undef, N)

    pos = SVector(zeros(T, 3)...)
    vy = SVector(0., 1., 0.)

    x = SAD*tan(ΔΘ)

    for i in eachindex(ϕg, n, p)
        src = SAD*SVector(sin(ϕg[i]), zero(T), cos(ϕg[i]))

        v = x*normalize(cross(pos-src, vy))

        pp = pos+v
        pm = pos-v

        pIc = closest_intersection(src, pos, mesh)
        pIp = closest_intersection(src, pp, mesh)
        pIm = closest_intersection(src, pm, mesh)

        p[i] = pIc
        n[i] = normalize(cross(pIp-pIm, vy))
    end

    LinearSurface(p, n)
end

function getplane(surf::LinearSurface, src::SVector{3})
    ϕg = atand(src[1], src[3])
    ϕg = mod(ϕg, 360)+1
    param = surf.params(ϕg)

    p = param[SVector(1, 2, 3)]
    n = param[SVector(4, 5, 6)]

    Plane(p, n)
end

function getSSD(surf::LinearSurface, pos, src)
    plane = getplane(surf, src)
    pI = intersection_point(plane, pos, src)
    pI === nothing && return Inf
    norm(pI-src)
end

function getdepth(surf::LinearSurface, pos, src)
    plane = getplane(surf, src)
    pI = intersection_point(plane, pos, src)
    pI === nothing && return NaN
    v = pI-pos
    sign(dot(plane.n, v))*norm(v)
end

#--- CylindricalSurface --------------------------------------------------------

"""
    CylindricalSurface

A planar external surface at a variable distance from the isocenter.
"""
struct CylindricalSurface{Ty<:AbstractVector, Tϕ<:AbstractVector, Tdist<:AbstractMatrix, TInterpolation} <: AbstractExternalSurface
    ϕ::Tϕ
    y::Ty
    distance::Tdist
    I::TInterpolation
    function CylindricalSurface(ϕ, y, rho)
        I = linear_interpolation((ϕ, y), rho)
        new{typeof(ϕ), typeof(y), typeof(rho), typeof(I)}(ϕ, y, rho, I)
    end
end

# Constructors

"""
    CylindricalSurface

Construct from a mesh.
"""
function CylindricalSurface(mesh::SimpleMesh; Δϕ°=2., Δy=2.)
    ϕ = (-180:Δϕ°:180)*pi/180

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
        pI = intersect_mesh(line, mesh)
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

    f(x) = distance_to_surface(x, surf, pos, src)

    sign(f(0)) == sign(f(1.)) && return Inf

    λ = find_zero(f, (0.5, 1.), AlefeldPotraShi()) #verbose=true
    λ*norm(src-pos)
end
getdepth(surf::CylindricalSurface, pos, src) = norm(pos - src) - getSSD(surf, pos, src)

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

function isinside(surf::CylindricalSurface, pos::AbstractVector{T}) where T<:Real
    x, y, z = pos

    (y<surf.y[1]||surf.y[end]<=y) && return false
    x^2+z^2 == zero(T) && return true

    ϕ = atan(x, z)
    x^2+z^2 < surf.I(ϕ, y)^2
end
