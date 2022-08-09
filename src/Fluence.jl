#
#   Fluence.jl
#
# Functions for creating fluence grids and computing fluence from beam-limiting
# devices.
#

export fluence, fluence!, bixel_grid, bixels_from_bld

export Bixel, position, width, area, subdivide

#--- Abstract Bixel -----------------------------------------------------------

"""
    AbstractBixel{T<:AbstractFloat}
"""
abstract type AbstractBixel{T<:AbstractFloat} end

"""
    Base.getindex(bixel::AbstractBixel, i::Int)

Get the centre position of the bixel.
"""
Base.getindex(bixel::AbstractBixel, i::Int) = bixel.position[i]

"""
    position(bixel::AbstractBixel)

Get the position of the bixel.
"""
position(bixel::AbstractBixel) = bixel.position

"""
    position(bixel::AbstractBixel, i::Int)

Get the ith coordinate of the position.
"""
position(bixel::AbstractBixel, i::Int) = bixel.position[i]

"""
    width(bixel::AbstractBixel)

Return the width of the bixel.
"""
width(bixel::AbstractBixel) = bixel.width

"""
    width(bixel::AbstractBixel, i::Int)

Return the width of the bixel along axis `i`.
"""
width(bixel::AbstractBixel, i::Int) = bixel.width[i]

"""
    area(bixel::AbstractBixel)

Return the area of the bixel.
"""
area(bixel::AbstractBixel) = prod(bixel.width)

"""
    subdivide(bixel::AbstractBixel, nx::Integer, ny::Integer)

Subdivide a bixel by specifing the number of partitions `nx` and `ny`.

Returns a grid of bixels.
"""
function subdivide(bixel::AbstractBixel, nx::Integer, ny::Integer)
    Δx, Δy = width(bixel)
    x, y = position(bixel)

    xsub = range(x-0.5*Δx, x+0.5*Δx, length=nx+1)
    xsub = 0.5*(xsub[1:end-1] .+ xsub[2:end])

    ysub = range(y-0.5*Δy, y+0.5*Δy, length=ny+1)
    ysub = 0.5*(ysub[1:end-1] .+ ysub[2:end])

    Bixel.(xsub, ysub', Δx/nx, Δy/ny)
end

"""
    subdivide(bixel::AbstractBixel{T}, δx::T, δy::T)

Subdivide by specifing widths `δx` and `δy`
"""
function subdivide(bixel::AbstractBixel{T}, δx::T, δy::T) where T<:AbstractFloat
    Δx, Δy = width(bixel)

    nx = ceil(Int, Δx/δx)
    ny = ceil(Int, Δy/δy)

    subdivide(bixel, nx, ny)
end

#--- Bixel --------------------------------------------------------------------

"""
    Bixel{T}

"""
struct Bixel{T} <: AbstractBixel{T}
    position::SVector{2, T}     # Centre position of the bixel
    width::SVector{2, T}   # Size of the bixel
end

Bixel(x::T, y::T, wx::T, wy::T) where T<:AbstractFloat = Bixel{T}(SVector(x, y), SVector(wx, wy))
Bixel(x::T, y::T, w::T) where T<:AbstractFloat = Bixel(x, y, w, w)
Bixel(x::T, w::T) where T<:AbstractFloat = Bixel(x, x, w, w)

function Bixel(x::AbstractVector{T}, w::AbstractVector{T}) where T<:AbstractFloat 
    Bixel{T}(x[1], x[2], w[1], w[2])
end

#--- Bixel Grid ---------------------------------------------------------------

"""
    snapped_range(x1, x2, Δ)

Create a range from x1 to x2 which is "snapped" to the step Δ.

Positions are "snapped" to the step value (e.g. a starting position of
x[1]-0.2Δx snaps to x[1]-Δx). The new range always includes the start and end
points of the original range

Examples:
- 0.1:1.:9.4 -> 0:1.:10.

"""
snapped_range(x1, x2, Δ) = Δ*(floor(Int, x1/Δ):ceil(Int, x2/Δ))

"""
    bixel_grid(x, y, Δx[, Δy])

Construct a grid of bixels.

Each axis starts at the first element (e.g. x[1]), runs to the last element
(x[end]), with uniform spacing (Δx). Same for y. Positions are "snapped" to the
spacing value (see `snapped_range` for details).

If Δy not specified, assumes Δy = Δx.
"""
function bixel_grid(x, y, Δx, Δy)
    x = snapped_range(x[1], x[end], Δx)
    y = snapped_range(y[1], y[end], Δy)

    Δx = @. x[2:end] - x[1:end-1]
    Δy = @. y[2:end] - y[1:end-1]

    x = @. 0.5*(x[1:end-1] + x[2:end])
    y = @. 0.5*(y[1:end-1] + y[2:end])

    Bixel.(SVector.(x, y'), SVector.(Δx, Δy'))
end

bixel_grid(x, y, Δ) = bixel_grid(x, y, Δ, Δ)

"""
    bixel_grid(jaws::Jaws, Δx[, Δy])

Uses the jaw positions to construct a bixel grid.

If Δy not specified, assumes Δy = Δx
"""
bixel_grid(jaws::Jaws, Δx, Δy) = bixel_grid(jaws.x, jaws.y, Δx, Δy)

bixel_grid(jaws::Jaws, Δ) = bixel_grid(jaws, Δ, Δ)

"""
    bixel_grid(x::AbstractRange, y::AbstractRange)

Uses the start and end positions and step of each range to construct the bixel grid.
"""
bixel_grid(x::AbstractRange, y::AbstractRange) = bixel_grid(x, y, step(x), step(y))

"""
    bixel_grid(mlc::MultiLeafCollimator, jaws::Jaws, Δx)

Grid that fits in an MLC and the jaws.

Bixel y widths are of the same width as the MLC leaf widths. Creates smaller widths
in the case where the jaws are halfway within a leaf width. Bixel x widths are
set by `Δx`.`
"""
function bixel_grid(mlc::MultiLeafCollimator, jaws::Jaws, Δx)

    iL, iU = subset_indices(mlc, jaws.y)

    y = Float64[]
    Δy = Float64[]

    for i = iL:iU
        yL, yU = mlc[i]

        yL = max(yL, jaws.y[1])
        yU = min(yU, jaws.y[2])

        push!(y, 0.5*(yL + yU))
        push!(Δy, yU - yL)

    end

    x = snapped_range(jaws.x[1], jaws.x[end], Δx)

    Bixel.(x, y', Δx, Δy')
end

"""
    bixels_from_bld(args::AbstractBeamLimitingDevice...)

Create bixels corresponding to the provided beam limiting devices.
""" bixels_from_bld

"""
    bixel_from_bld(jaws::Jaws)

From Jaws.
"""
function bixels_from_bld(jaws::Jaws{T}) where T<:AbstractFloat

    x = T(0.5)*(jaws.x[1]+jaws.x[2])
    y = T(0.5)*(jaws.y[1]+jaws.y[2])

    Δx = jaws.x[2] - jaws.x[1]
    Δy = jaws.y[2] - jaws.y[1]

    [Bixel(x, y, Δx, Δy)]
end

"""
    bixels_from_bld(mlcx, mlc::MultiLeafCollimator, jaws::Jaws)

From MultiLeafCollimator and Jaws.
"""
function bixels_from_bld(mlcx::AbstractMatrix{T}, mlc::MultiLeafCollimator, jaws::Jaws{T}) where T<:AbstractFloat

    bixels = Bixel{T}[]

    @inbounds for i in eachindex(mlc)
        xL, xU = mlcx[:, i]
        yL, yU = mlc[i]

        xL = max(xL, jaws.x[1])
        xU = min(xU, jaws.x[2])

        yL = max(yL, jaws.y[1])
        yU = min(yU, jaws.y[2])

        Δy = yU - yL
        Δx = xU - xL

        if(Δy > zero(T) && Δx > zero(T))
            x = T(0.5)*(xL + xU)
            y = T(0.5)*(yL + yU)

            push!(bixels, Bixel(x, y, Δx, Δy))
        end
    end
    bixels

end

#--- Fluence in a single bixel -----------------------------------------------------------------------------------------

"""
    overlap(x, Δx, xB, xA)

Compute the overlapping length of a bixel spanning `x`-`w`/2->`x`+`w`/2 and the
length between `xL` and `xU`, normalised to the length of the bixel. If the
Bixel fully within range (xL<=x-w/2 && x+w/2<=xU), return 1. If the bixel is
fully outside the range (xU<=x-w/2 || x+w/2<=xL), return 0.

"""
function overlap(x::T, w::T, xL::T, xU::T) where T<:AbstractFloat
    hw = T(0.5)*w
    max(zero(T), min(x + hw, xU) - max(x - hw, xL))/w
end

"""
    fluence_from_rectangle(bixel::Bixel, xlim, ylim)

Compute the fluence of a rectangle with edges at `xlim` and `ylim` on a bixel.
"""
function fluence_from_rectangle(bixel::Bixel, xlim, ylim)
    overlap(position(bixel, 1), width(bixel, 1), xlim[1], xlim[2])*overlap(position(bixel, 2), width(bixel, 2), ylim[1], ylim[2])
end

#--- Computing Fluence -------------------------------------------------------------------------------------------------

"""
    fluence(bixel::AbstractBixel, bld::AbstractBeamLimitingDevice, args...)

Compute the fluence of `bixel` from beam limiting device (e.g. an MLC or jaws).
"""
function fluence(bixel::AbstractBixel, bld::AbstractBeamLimitingDevice, args...) end

"""
    fluence(bixels::AbstractArray{<:AbstractBixel}, bld::AbstractBeamLimitingDevice, args...)

Compute the fluence on a collection of bixels.

Broadcasts over the specific `fluence(bixel, ...)` method for the provided beam
limiting device.
*e.g.*: `fluence(bixels, jaws)`, `fluence(bixels, mlcx, mlc)`
"""
fluence(bixels::AbstractArray{<:AbstractBixel}, args...) = fluence.(bixels, Ref.(args)...)

fluence!(Ψ::AbstractArray{<:AbstractFloat}, bixels::AbstractArray{<:AbstractBixel}, args...) = Ψ .= fluence.(bixels, Ref.(args)...)


"""
    fluence(bixels::AbstractArray{<:AbstractBixel}, index::AbstractArray{Int}, args...)

Allows precomputation of the location of the bixel in relation to the beam limiting device.

In the case of an MLC, the index is the leaf index that contains that bixel.

Requires `fluence(bixel::AbstractBixel, index::Int, args...)` to be defined for
the particular beam limiting device
"""
fluence(bixels::AbstractArray{<:AbstractBixel}, index::AbstractArray{Int}, args...) = fluence.(bixels, index, Ref.(args)...)

function fluence!(Ψ::AbstractArray{<:AbstractFloat}, bixels::AbstractArray{<:AbstractBixel}, index::AbstractArray{Int}, args...)
    Ψ .= fluence.(bixels, index, Ref.(args)...)
end



#--- Computing Fluence from Jaws ---------------------------------------------------------------------------------------

"""
    fluence(bixel::Bixel, jaws::Jaws)

From the Jaws.
"""
fluence(bixel::Bixel, jaws::Jaws) = fluence_from_rectangle(bixel, jaws.x, jaws.y)

#--- Fluence from an MLC Aperture --------------------------------------------------------------------------------------

"""
    fluence(bixel::Bixel, mlcx, mlc::MultiLeafCollimator)

From an MLC aperture.
"""
function fluence(bixel::AbstractBixel{T}, mlcx, mlc::MultiLeafCollimator) where T<:AbstractFloat

    hw = 0.5*width(bixel, 2)
    i1 = max(1, locate(mlc, bixel[2]-hw))
    i2 = min(length(mlc), locate(mlc, bixel[2]-hw))

    Ψ = zero(T)
    @inbounds for j=i1:i2
        Ψ += fluence_from_rectangle(bixel, (@view mlcx[:, j]), (@view mlc[j]))
    end
    Ψ
end

"""
    fluence(bixel::Bixel, index::Int, mlcx)

From an MLC aperture using a given leaf index.

This method assumes the bixel is entirely within the `i`th leaf track, and does
not overlap with other leaves. Does not check whether these assumptions are true.
"""
function fluence(bixel::AbstractBixel{T}, index::Int, mlcx) where T<:AbstractFloat
    overlap(position(bixel, 1), width(bixel, 1), mlcx[1, index], mlcx[2, index])
end

#--- Moving Aperture fluences ------------------------------------------------------------------------------------------

"""
    fluence(bixel::AbstractBixel{T}, mlcx1, mlcx2, mlc::MultiLeafCollimator)

From MLC leaf positions which move from `mlcx1` to `mlcx2`.

Computes the time-weighted fluence as the MLC moves from position `mlcx1` to `mlcx2`.
Assumes the MLC leaves move in a straight line.
"""
function fluence(bixel::AbstractBixel{T}, mlcx1, mlcx2, mlc::MultiLeafCollimator) where T<:AbstractFloat
    j = max(1, locate(mlc, bixel[2]))
    xL = bixel[1] - 0.5*width(bixel, 1)
    xU = bixel[1] + 0.5*width(bixel, 1)

    ΨB = fluence_onesided(mlcx1[1, j], mlcx2[1, j], xL, xU)
    ΨA = fluence_onesided(mlcx1[2, j], mlcx2[2, j], xL, xU)

    max(zero(T), ΨB - ΨA) # This is needed for cases where xA < xB
end

"""
    fluence_onesided(xs, xf, xL, xU)

Compute fluence for 1 leaf position, assuming the other is infinitely far away

Computes the fluence for a leaf trajectory from `xs` to `xf`, over a bixel from
`xL` to `xU`. The aperture is considered open to the right of the leaf position
(*i.e.* leaf on the B bank). Assumes the bixel is fully in the leaf track.
"""
function fluence_onesided(xs, xf, xL, xU)

    # If leaf positions are the same, compute static fluence
    xs == xf && return (xU - minmax(xs, xL, xU))/(xU-xL)


    # Otherwise compute moving fluence

    x1, x2 = min(xs, xf), max(xs, xf) # If xf < xs, swap

    xl = max(x1, xL) # Get left, centre and right segment positions
    xc = min(x2, xU)
    xr = xU

    tl = leaf_trajectory(xl, x1, x2) # Compute the height of each segment pos.
    tr = leaf_trajectory(xr, x1, x2)
    
    A1 = 0.5*(tl + tr)*(xc-xl)  # Calculate trapezium segment are
    A2 = xr-xc  # Calculate the remaining rectangular area
    (A1 + A2)/(xU-xL)   # Add and scale

end

"""
    minmax(x, l, u)

Return `x` if `l<=x<=u`, `l` if `x<l` or `u` if `u<x`
"""
minmax(x, l, u) = max(l, min(x, u)) # Ensures l<=x<=u

"""
    leaf_trajectory(x, x1, x2)

Compute the height of position `x` between `(x1, 0)` and `(x2, 1)`

Used in `intersection_area`.
"""
leaf_trajectory(x, x1, x2) = minmax((x-x1)/(x2-x1), 0, 1)
