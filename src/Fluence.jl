#
#   Fluence.jl
#
# Functions for creating fluence grids and computing fluence from beam-limiting
# devices.
#

export fluence, fluence!, bixel_grid, bixels_from_bld

export Bixel, getcenter, getwidth, getedge, getarea, subdivide

#--- Abstract Fluence Element -------------------------------------------------

"""
    AbstractFluenceElement
"""
abstract type AbstractFluenceElement end

#--- Abstract Bixel -----------------------------------------------------------

"""
    AbstractBixel
"""
abstract type AbstractBixel <: AbstractFluenceElement end

#--- Bixel --------------------------------------------------------------------

"""
    Bixel{T}

"""
struct Bixel{T} <: AbstractBixel
    position::SVector{2, T}     # Centre position of the bixel
    width::SVector{2, T}   # Size of the bixel
end

Bixel(x::T, y::T, wx::T, wy::T) where T<:AbstractFloat = Bixel{T}(SVector(x, y), SVector(wx, wy))
Bixel(x::T, y::T, w::T) where T<:AbstractFloat = Bixel(x, y, w, w)
function Bixel(x::AbstractVector, w::AbstractVector)
    x, w = promote(x, w)
    Bixel(SVector(x...), SVector(w...))
end

# function Bixel(x::AbstractVector{T}, w::AbstractVector{T}) where T<:AbstractFloat 
#     Bixel{T}(x[1], x[2], w[1], w[2])
# end


"""
    Base.getindex(bixel::Bixel, i::Int)

Get the centre position of the bixel.
"""
Base.getindex(bixel::Bixel, i::Int) = bixel.position[i]

"""
    getcenter(bixel::Bixel)

Get the position of the bixel.
"""
getcenter(bixel::Bixel) = bixel.position

"""
    getcenter(bixel::Bixel, i::Int)

Get the ith coordinate of the position.
"""
getcenter(bixel::Bixel, i::Int) = bixel.position[i]

"""
    getwidth(bixel::Bixel)

Return the width of the bixel.
"""
getwidth(bixel::Bixel) = bixel.width

"""
    getwidth(bixel::Bixel, i::Int)

Return the width of the bixel along axis `i`.
"""
getwidth(bixel::Bixel, i::Int) = bixel.width[i]

"""
    getedge(bixel::Bixel[, dim::Int])

Return the lower edge of the bixel.

Can specify a `dim` for which dimension (x or y).
"""
getedge(bixel::Bixel, args...) = getcenter(bixel, args...) - 0.5*getwidth(bixel, args...)

"""
    getarea(bixel::Bixel)

Return the area of the bixel.
"""
getarea(bixel::Bixel) = prod(bixel.width)

"""
    subdivide(bixel::Bixel, nx::Integer, ny::Integer)

Subdivide a bixel by specifing the number of partitions `nx` and `ny`.

Returns a grid of bixels.
"""
function subdivide(bixel::Bixel, nx::Integer, ny::Integer)
    Δx, Δy = getwidth(bixel)
    x, y = getcenter(bixel)

    xsub = range(x-0.5*Δx, x+0.5*Δx, length=nx+1)
    xsub = 0.5*(xsub[1:end-1] .+ xsub[2:end])

    ysub = range(y-0.5*Δy, y+0.5*Δy, length=ny+1)
    ysub = 0.5*(ysub[1:end-1] .+ ysub[2:end])

    Bixel.(xsub, ysub', Δx/nx, Δy/ny)
end

"""
    subdivide(bixel::Bixel{T}, δx::T, δy::T)

Subdivide by specifing widths `δx` and `δy`
"""
function subdivide(bixel::Bixel{T}, δx::T, δy::T) where T<:AbstractFloat
    Δx, Δy = getwidth(bixel)

    nx = ceil(Int, Δx/δx)
    ny = ceil(Int, Δy/δy)

    subdivide(bixel, nx, ny)
end

#--- Bixel Grid ---------------------------------------------------------------

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
    bixel_grid(mlc::AbstractMultiLeafCollimator, jaws::Jaws, Δx)

Grid that fits in an MLC and the jaws.

Bixel y widths are of the same width as the MLC leaf widths. Creates smaller widths
in the case where the jaws are halfway within a leaf width. Bixel x widths are
set by `Δx`.`
"""
function bixel_grid(mlc::AbstractMultiLeafCollimator, jaws::Jaws, Δx)

    iL = locate(mlc, jaws.y[1])
    iU = locate(mlc, jaws.y[2])

    y = getedges(mlc)[iL:iU]
    y = minmax.(y, jaws.y[1], jaws.y[2])
    Δy = diff(y)

    x = snapped_range(jaws.x[1], jaws.x[end], Δx)
    x = minmax.(x, jaws.x[1], jaws.x[2])
    Δx = diff(x)

    x = 0.5*(x[1:end-1] + x[2:end])
    y = 0.5*(y[1:end-1] + y[2:end])

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
function bixels_from_bld(mlc::AbstractMultiLeafCollimator, jaws::Jaws{T}; Δx=5., Δy=5.) where T<:AbstractFloat

    bixels = Bixel{T}[]

    @inbounds for i in eachindex(mlc)
        ((xL, xU), (yL, yU)) = mlc[i]
        
        xL = max(xL, jaws.x[1])
        xU = min(xU, jaws.x[2])

        yL = max(yL, jaws.y[1])
        yU = min(yU, jaws.y[2])

        if xU - xL > zero(T)
            x = snapped_range(xL, xU, Δx)
            # x = max.(xL, min.(x, xU))
            xc = 0.5*(x[2:end] + x[1:end-1])
            δx = x[2:end] - x[1:end-1]

            y = 0.5*(yL + yU)

            b = Bixel.(xc, y, δx, Δy)

            append!(bixels, vec(b))
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
function overlap(x::Number, w::Number, xL::Number, xU::Number)
    hw = 0.5*w
    max(0., min(x + hw, xU) - max(x - hw, xL))/w
end

"""
    fluence_from_rectangle(bixel::Bixel, xlim, ylim)

Compute the fluence of a rectangle with edges at `xlim` and `ylim` on a bixel.
"""
function fluence_from_rectangle(bixel::Bixel, xlim, ylim)
    overlap(getcenter(bixel, 1), getwidth(bixel, 1), xlim[1], xlim[2])*overlap(getcenter(bixel, 2), getwidth(bixel, 2), ylim[1], ylim[2])
end

#--- Computing Fluence -------------------------------------------------------------------------------------------------

"""
    fluence(bixel::AbstractBixel, bld::AbstractBeamLimitingDevice, args...)

Compute the fluence of `bixel` from beam limiting device (e.g. an MLC or jaws).
""" fluence(bixel::AbstractBixel, bld::AbstractBeamLimitingDevice, args...)

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
    fluence(bixel::Bixel, mlc::MultiLeafCollimator)

From an MLC aperture.
"""
function fluence(bixel::Bixel{T}, mlc::MultiLeafCollimator) where T<:AbstractFloat

    hw = 0.5*getwidth(bixel, 2)
    i1 = max(1, locate(mlc, bixel[2]-hw))
    i2 = min(length(mlc), locate(mlc, bixel[2]-hw))

    Ψ = zero(T)
    @inbounds for j=i1:i2
        Ψ += fluence_from_rectangle(bixel, mlc[j]...)
    end
    Ψ
end

"""
    fluence(bixel::Bixel, index::Int, mlcx)

From an MLC aperture using a given leaf index.

This method assumes the bixel is entirely within the `i`th leaf track, and does
not overlap with other leaves. Does not check whether these assumptions are true.
"""
function fluence(bixel::AbstractBixel, index::Int, mlcx::AbstractMatrix)
    overlap(getcenter(bixel, 1), getwidth(bixel, 1), mlcx[1, index], mlcx[2, index])
end

#--- Moving Aperture fluences ------------------------------------------------------------------------------------------

"""
    fluence(bixel::Bixel{T}, mlc1::MultiLeafCollimator, mlc2::MultiLeafCollimator)

From an MLC aperture sequence.
"""
function fluence(bixel::Bixel{T}, mlc1::MultiLeafCollimator, mlc2::MultiLeafCollimator) where T<:AbstractFloat

    hw = 0.5*getwidth(bixel, 2)
    i1 = max(1, locate(mlc1, bixel[2]-hw))
    i2 = min(length(mlc1), locate(mlc1, bixel[2]-hw))

    Ψ = zero(T)
    @inbounds for j=i1:i2
        Ψ += fluence_from_moving_aperture(bixel, getpositions(mlc1, j), getpositions(mlc2, j))
    end
    Ψ
end

"""
    fluence(bixel::Bixel{T}, index::Int, mlcx1, mlcx2)

From an MLC aperture sequence using a given leaf index.
"""
function fluence(bixel::AbstractBixel, index::Int, mlcx1::AbstractMatrix, mlcx2::AbstractMatrix)
    fluence_from_moving_aperture(bixel, (@view mlcx1[:, index]), (@view mlcx2[:, index]))
end

"""
    fluence_from_moving_aperture(bixel::Bixel{T}, mlcx1, mlcx2)

From MLC leaf positions which move from `mlcx1` to `mlcx2`.

Computes the time-weighted fluence as the MLC moves from position `mlcx1` to `mlcx2`.
Assumes the MLC leaves move in a straight line.
"""
function fluence_from_moving_aperture(bixel::Bixel{T}, mlcx1, mlcx2) where T<:AbstractFloat
    xL = bixel[1] - 0.5*getwidth(bixel, 1)
    xU = bixel[1] + 0.5*getwidth(bixel, 1)

    ΨB = fluence_onesided(mlcx1[1], mlcx2[1], xL, xU)
    ΨA = fluence_onesided(mlcx1[2], mlcx2[2], xL, xU)

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
    leaf_trajectory(x, x1, x2)

Compute the height of position `x` between `(x1, 0)` and `(x2, 1)`

Used in `intersection_area`.
"""
leaf_trajectory(x, x1, x2) = minmax((x-x1)/(x2-x1), 0, 1)
