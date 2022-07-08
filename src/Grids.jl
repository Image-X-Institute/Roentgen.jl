
import Base.+, Base.-, Base.size, Base.getindex, Base.length, Base.step, Base.show, Base.eachindex, Base.lastindex

export GridUniform, GridUniform2D, locate, spacing, interp, getx, getxc, gety, getyc

abstract type AbstractGrid{N} <: AbstractArray{AbstractFloat,N} end # where {T<:AbstractFloat, N<:Integer} end


## 1D Uniform Grid

struct GridUniform{T} <: AbstractGrid{1}
    range::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
end

# Basic Methods
Base.length(grid::GridUniform) = length(grid.range)-1
Base.size(grid::GridUniform) = (length(grid),)

Base.getindex(grid::GridUniform, i::Integer) = getx(grid, i)


Base.iterate(grid::GridUniform) = iterate(grid.range)
Base.iterate(grid::GridUniform, state) = iterate(grid.range, state)

locate(grid::GridUniform{T}, xi::T) where T<:AbstractFloat = locate(getx(grid), xi)
spacing(grid::GridUniform) = step(grid.range)

getx(grid::GridUniform) = grid.range
getx(grid::GridUniform, i::Integer) = grid.range[i]

getxc(grid::GridUniform) = grid.range[1:end-1] .+ 0.5*spacing(grid)
getxc(grid::GridUniform, i::Integer) = getx(grid, i) + 0.5*spacing(grid)


Base.:+(grid::GridUniform, x::Number) = GridUniform(getx(grid) .+ x)
Base.:-(grid::GridUniform, x::Number) = GridUniform(getx(grid) .- x)

Base.show(io::IO, grid::GridUniform) = show(io, grid.range)

# 2D Uniform Grid

struct GridUniform2D{T} <: AbstractGrid{2}
    x::GridUniform{T}
    y::GridUniform{T}
end

Base.show(io::IO, grid::GridUniform2D) = show(io, (grid.x, grid.y))

Base.length(grid::GridUniform2D) = length(grid.x)*length(grid.y)
Base.size(grid::GridUniform2D) = (size(grid.x)..., size(grid.y)...)

Base.eachindex(grid::GridUniform2D) = Base.OneTo(length(grid))

Base.getindex(grid::GridUniform2D, i::Integer, j::Integer) = getx(grid, i), gety(grid, j)
Base.getindex(grid::GridUniform2D, index::Integer) = grid[CartesianIndices(size(grid))[index]]
Base.getindex(grid::GridUniform2D, index::CartesianIndex) = grid[index[1], index[2]]

Base.:+(grid::GridUniform2D, x) = GridUniform2D(grid.x + x[1], grid.y + x[2])
Base.:-(grid::GridUniform2D, x) = GridUniform2D(grid.x - x[1], grid.y - x[2])

function GridUniform2D(x::T, y::T) where T<:StepRangeLen
    gridx = GridUniform(x)
    gridy = GridUniform(y)
    GridUniform2D(gridx, gridy)
end

function GridUniform2D(x)
    grid = GridUniform(x)
    GridUniform2D(grid, grid)
end

locate(grid::GridUniform2D, x, y) = CartesianIndex(locate(grid.x, x), locate(grid.y, y))
locate(grid::GridUniform2D, p) = locate(grid, p[1], p[2])

spacing(grid::GridUniform2D) = spacing(grid.x), spacing(grid.y)

getx(grid::GridUniform2D) = getx(grid.x)
getx(grid::GridUniform2D, i::Integer) = getx(grid.x, i)

gety(grid::GridUniform2D) = getx(grid.y)
gety(grid::GridUniform2D, i::Integer) = getx(grid.y, i)

getxc(grid::GridUniform2D) = getxc(grid.x)
getxc(grid::GridUniform2D, i::Integer) = getxc(grid.x, i)

getyc(grid::GridUniform2D) = getxc(grid.y)
getyc(grid::GridUniform2D, i::Integer) = getxc(grid.y, i)

## Interpolation

### Linear Interpolation
interp(grid::GridUniform, f::AbstractVector, xi) = interp(getxc(grid), f, xi)

### Bilinear Interpolation
function interp(grid::GridUniform2D, f, xi, yi)
    interp(grid.x.range[1:end-1], grid.y.range[1:end-1], f, xi, yi)
end
