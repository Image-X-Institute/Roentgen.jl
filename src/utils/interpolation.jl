#
#   Linear Interpolation Functions
#

export locate, interp

#--- Utility Functions --------------------------------------------------------

"""
    scale_to_cell(x1, x2, xi)

Scale the position `xi` within the positions `x1` and `x2`.
"""
scale_to_cell(x1, x2, xi) = (xi - x1)/(x2 - x1)

"""
    minmax(xmin, xmax, xi)

Limit the position `xi` between `xmin` and `xmax`.
"""
minmax(xmin, xmax, x) = max(xmin, min(xmax, x))

#--- Non-Uniform Grid Location ------------------------------------------------

"""
    locate(xg::AbstractVector{T}, xi::T) where T<:Real

Locates the grid index for position `x` within a grid `xg`. Does not limit out
of bounds indices.
"""
locate(xg::AbstractVector{T}, xi::T) where T<:Real = searchsortedlast(xg, xi)

#--- Uniform Grid Location ----------------------------------------------------

"""
    locate(xg::AbstractRange{T}, xi::T) where T<:Real

When a range is provided, uniform spacing can be assumed.
"""
locate(xg::AbstractRange{T}, xi::T) where T<:Real = locate(xg[1], step(xg), xi)

"""
    locate(x::T, start::T, step::T) where T<:AbstractFloat

Locate the grid index for position `xi`

Specify the grid starting position (`start`) and grid spacing (`step`).
"""
@inline locate(start, step, xi) = floor(Int, (xi - start)/step) + 1

#--- 2D Location --------------------------------------------------------------

locate(xg, yg, xi, yi) = CartesianIndex(locate(xg, xi), locate(yg, yi))

#--- Grid Interpolation -------------------------------------------------------

"""
    interp(xg::AbstractVector{T}, fg::AbstractVector{T}, xi) where T<:Real

Interpolate `xi` within uniform grid positions `xg` and grid values `fg`
"""
Base.@propagate_inbounds function interp(xg::AbstractVector{T}, fg::AbstractVector{T}, xi) where T<:Real
    LinearInterpolation(xg, fg, extrapolation_bc = Flat())(xi)
end

"""
    interp(f1, f2, α)

Interpolate using given values `f1` and `f2` and the scaled distance from them `α`
"""
interp(f1, f2, α) = (1 - α)*f1 + α*f2

#--- Bilinear Interpolation ---------------------------------------------------

"""
    interp(xg::AbstractVector, yg::AbstractVector, fg::AbstractMatrix, xi, yi)

Bilinear interpolation at position `xi,yi`, on grid `xg-yg` with values `fg`
"""
Base.@propagate_inbounds function interp(xg::AbstractVector, yg::AbstractVector, fg::AbstractMatrix, xi, yi)
    LinearInterpolation((xg, yg), fg, extrapolation_bc = Flat())(xi, yi)
end

