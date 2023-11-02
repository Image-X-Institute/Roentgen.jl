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

getdepth(surf::ConstantSurface, pos, src) = norm(pos-src)-getSSD(surf, pos, src)

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

getdepth(surf::PlaneSurface, pos, src) = norm(pos-src)-getSSD(surf, pos, src)

#--- MeshSurface ---------------------------------------------------------------

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
    function MeshSurface(mesh::Partition, boxes)
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

"""
    LinearSurface(I::AbstractInterpolation})

Linear approximation of a general surface by gantry angle

Constructed with an interpolator which returns the plane position and normal
vectors at a given gantry angle.
"""
struct LinearSurface{T<:AbstractInterpolation} <: AbstractExternalSurface
    params::T

    function LinearSurface(params::AbstractVector)
        @assert length(params)==361 "Distance must be supplied at every degree"
        I = interpolate(params, BSpline(Linear()))
        new{typeof(I)}(I)
    end
    LinearSurface(I::AbstractInterpolation) = new{typeof(I)}(I)
end

"""
    LinearSurface(params::AbstractVector)

Constructed with a vector of 6 element parameters corresponding to the plane 
position and normal at gantry angles.
"""
LinearSurface

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
function LinearSurface(mesh::SimpleMesh{3, T}; SAD=T(1000.), Δϕ=deg2rad(1)) where {T<:Real}
    N = 361
    ϕg = 2π*range(0, 1, length=N)
    n = Vector{SVector{3, T}}(undef, N)
    p = Vector{SVector{3, T}}(undef, N)

    pos = SVector(zeros(T, 3)...)
    vy = SVector(0., 1., 0.)

    x = SAD*tan(Δϕ)

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

_assert_cylsurf_size(ϕ, y, rho) = @assert size(rho) == (length(ϕ), length(y))
function _warn_cylsurf_rho(rho)
    if any(@. !isfinite(rho) || rho<0)
        @warn "rho contains negative values, Inf or NaNs. Ensure that the surface
        is well defined in cylindrical coordinate system with the axial direction
        along the y axis. This may be fixed by lowering the resolution in the `y`
        and `ϕ` axes"
    end
end

"""
    CylindricalSurface(y::TRange, ϕ::TRange, rho::TInterp)

Surface stored on a cylindrical-polar grid.
"""
struct CylindricalSurface{TVector<:AbstractVector, TRange<:AbstractRange, TInterp<:AbstractInterpolation} <: AbstractExternalSurface
    ϕ::TRange
    y::TRange
    rho::TInterp
    center::TVector
    function CylindricalSurface(ϕ::TRange, y::TRange, rho::TInterp, center::TVector) where {
        TVector<:AbstractVector, TRange<:AbstractRange, TInterp<:AbstractInterpolation}
        _assert_cylsurf_size(ϕ, y, rho)
        _warn_cylsurf_rho(rho)
        new{TVector, TRange, TInterp}(ϕ, y, rho, center)
    end
end

# Constructors


"""
    CylindricalSurface(ϕ::AbstractVector, y::AbstractVector, rho::AbstractMatrix)

Constructed using vectors for ϕ, y, and rho.
"""
function CylindricalSurface(ϕ::AbstractVector, y::AbstractVector, rho::AbstractMatrix, center::AbstractVector)

    _assert_cylsurf_size(ϕ, y, rho)

    ϕI = 0:minimum(diff(ϕ)):2π
    yI = y[1]:minimum(diff(y)):y[end]

    rhoI = linear_interpolation((ϕ, y), rho).(ϕI, yI')

    I = interpolate(rhoI, BSpline(Linear()))
    CylindricalSurface(ϕI, yI, I, center)
end

"""
    CylindricalSurface(mesh::SimpleMesh, y::AbstractRange[, nϕ=181])

Construct from `mesh` over axial axis `y`.

Defaults to a 2° azimuthal spacing (nϕ=181).
"""
function CylindricalSurface(mesh::SimpleMesh, y::AbstractRange, nϕ::Int=181)
    # Compute azimuthal range
    ϕ = range(0., 2π, length=nϕ)

    # Set the length of the line used in computing the intersect
    L = 2*diagonal(boundingbox(mesh))

    # Center the y axis on the mesh center
    center = coordinates(centroid(mesh))
    y = y .- center[2]

    rho = fill(Inf, length(ϕ), length(y))

    for j in eachindex(y), i in eachindex(ϕ[1:end-1])

        # Create line along the axial axis, spanning radially out
        p1 = SVector(0., y[j], 0.) + center
        p2 = SVector(L*sin(ϕ[i]), y[j], L*cos(ϕ[i])) + center

        # Find the closest intersection along this line on the mesh
        pI = closest_intersection(p1, p2, mesh)

        # Move back into mesh coords
        x, _, z = pI-center

        # Compute distance from axis (radius)
        rho[i, j] = √(x^2+z^2)

    end
    # ϕ=0° and ϕ=360° are the same point
    rho[end, :] .= rho[1, :]

    # Interpolate, and return CylindricalSurface
    I = interpolate(rho, BSpline(Linear()))
    CylindricalSurface(ϕ, y, I, center)
end

"""
    CylindricalSurface(mesh::SimpleMesh, Δy::Real[, nϕ=181])

Construct from `mesh` over with axial spacing `y`.

Uses the mesh bounds to compute the axial range.
Defaults to a 2° azimuthal spacing (nϕ=181).
"""
function CylindricalSurface(mesh::SimpleMesh, Δy::Real, args...; kwargs...)
    box = boundingbox(mesh)
    y₀ = coordinates(minimum(box))[2]
    y₁ = coordinates(maximum(box))[2]
    y = snapped_range(y₀, y₁, Δy)[2:end-1]

    CylindricalSurface(mesh, y, args...; kwargs...)
end

# Methods

function _interp(surf::CylindricalSurface, ϕ, y)
    ϕg = surf.ϕ
    i = mod2pi(ϕ)/step(ϕg)+1

    yg = surf.y
    j = (y-first(yg))/step(yg)
    j = clamp(j+1, 1, length(yg))

    surf.rho(i, j)
end

function _interp(surf::CylindricalSurface, r)
    ϕ = atan(r[1], r[3])
    _interp(surf, ϕ, r[2])
end

function _distance_to_surface(λ, surf, pos, src)
    r = pos + λ*(src-pos) - surf.center

    rho = _interp(surf, r)

    idx = SVector(1, 3)
    rho^2 - dot(r[idx], r[idx])
end

function getdepth(surf::CylindricalSurface, pos::AbstractVector{T}, src::AbstractVector{T}) where T<:Real

    f(x) = _distance_to_surface(x, surf, pos, src)

    lims = (zero(T), one(T))

    sign(f(lims[1])) == sign(f(lims[2])) && return Inf

    λ = find_zero(f, lims, AlefeldPotraShi())
    λ*norm(src-pos)
end

function extent(surf::CylindricalSurface)
    ϕ = surf.ϕ
    rho = Interpolations.coefficients(surf.rho)
    xc, yc, zc = surf.center

    x = @. xc + rho*sin(ϕ)
    y = @. yc + surf.y
    z = @. zc + rho*cos(ϕ)
    SVector(minimum(x), minimum(y), minimum(z)), SVector(maximum(x), maximum(y), maximum(z))
end

function isinside(surf::CylindricalSurface, pos::AbstractVector{T}) where T<:Real
    x, y, z = pos-surf.center

    (y<surf.y[1]||surf.y[end]<=y) && return false

    ρ² = x^2+z^2
    ρ² == zero(T) && return true

    ρ²<_interp(surf, pos-surf.center)^2
end
