import Base.+, Base.-, Base.getindex, Base.lastindex, Base.size, Base.length
import Base.(==), Base.view
using StaticArrays

export Jaws, MultiLeafCollimator, getx, gety, inaperture, centerposition, edgeposition, extract_subset

"""
    Abstract type for Beam Limiting Devices are based on

`AbstractBeamLimitingDevice` contains the types for MultiLeafCollimators, Jaws, etc.
"""
abstract type AbstractBeamLimitingDevice end

abstract type AbstractMultiLeafCollimator end

#--- MultiLeafCollimator ------------------------------------------------------

"""
    MultiLeafCollimator

`MultiLeafCollimator` stores the leaf positions and edges of an MLC.

The 2 by `n` leaf positions are stored in the `positions` matrix. The leaf position
boundaries are stored in the `edges` vector, which is 1 element longer than the
number of leaf tracks, `n`.

Indexing a `MultiLeafCollimator` will either return the x and y bounds of the leaf
track, or another `MultiLeafCollimator`. It also supports the `@view` macro.
Iteration goes through each leaf track, returning the lower and upper boundaries.

A `MultiLeafCollimator` can also be shifted through addition or subtraction, and
scaled with multplication and division.

The `locate` method is implemented to return the leaf index containing edge position.
"""
struct MultiLeafCollimator{Tpos<:AbstractMatrix, Tedge<:AbstractVector} <: AbstractBeamLimitingDevice
    positions::Tpos
    edges::Tedge
    n::Int
    function MultiLeafCollimator(positions, edges)
        n = length(edges)-1
        @assert size(positions, 1) == 2 "Number of leaf positions per track != 2"
        @assert size(positions, 2) == n "Length of positions and edges do not match"
        new{typeof(positions), typeof(edges)}(positions, edges, n)
    end
end

# Constructors

"""
    MultiLeafCollimator(n::Int, Δy::Real)

Construct an MLC with `n` number of leaves and of leaf width `Δy`, centered
about zero.
"""
MultiLeafCollimator(n::Int, Δy::Real) = MultiLeafCollimator(zeros(2, n), Δy*(-0.5*n:0.5*n))

"""
    MultiLeafCollimator(n::Int, Δy::Real)

Construct an MLC with leaf edges.
"""
function MultiLeafCollimator(edges)
    n = length(edges)-1
    MultiLeafCollimator(zeros(2, n), y)
end

# Size and Length
Base.length(mlc::MultiLeafCollimator) = mlc.n
Base.size(mlc::MultiLeafCollimator) = (length(mlc),)

Base.firstindex(mlc::MultiLeafCollimator) = 1
Base.lastindex(mlc::MultiLeafCollimator) = mlc.n
Base.eachindex(mlc::MultiLeafCollimator) = Base.OneTo(mlc.n)

function Base.show(io::IO, mlc::MultiLeafCollimator)
    maxdigits = 6
    print(size(mlc.positions, 1), "x", size(mlc.positions, 2), " MultiLeafCollimator")
    for j =1:2
        println(io)
        for i in eachindex(mlc)

            msg = string(round(mlc.positions[j, i]; digits=maxdigits))

            colwidth = maximum(@. length(string(round(mlc.positions[:, i]; digits=maxdigits))))
            msgwidth = length(msg)

            print(io, " ", msg, repeat(" ", colwidth-msgwidth+1))
        end
    end
end

# Indexing
Base.getindex(mlc::MultiLeafCollimator, i::Int) = mlc.positions[:, i], mlc.edges[i:i+1]
Base.getindex(mlc::MultiLeafCollimator, i::UnitRange{Int}) = MultiLeafCollimator(mlc.positions[:, i], mlc.edges[i[1]:i[end]+1])

Base.view(mlc::MultiLeafCollimator, i::Int) = (@view mlc.positions[:, i]), (@view mlc.edges[i:i+1])
Base.view(mlc::MultiLeafCollimator, i::UnitRange{Int}) = MultiLeafCollimator((@view mlc.positions[:, i]), (@view mlc.edges[i[1]:i[end]+1]))

# Iteration
Base.iterate(mlc::MultiLeafCollimator) = (mlc[1], 1)
Base.iterate(mlc::MultiLeafCollimator, state) = state>=length(mlc) ? nothing : (mlc[state+1], state+1)


# Required for tests.
function Base.:(==)(mlc1::MultiLeafCollimator, mlc2::MultiLeafCollimator)
    mlc1.positions == mlc2.positions && mlc1.edges == mlc2.edges && mlc1.n == mlc2.n
end

#= Shifting the MLC

    Adds ability to move and scale the MLC
    From https://docs.julialang.org/en/v1/manual/metaprogramming/#Code-Generation
=#
for op in (:+, :-, :*, :/)
    eval(quote
        function Base.$op(mlc::MultiLeafCollimator, x::Union{AbstractVector,Tuple})
            MultiLeafCollimator( ($op).(mlc.positions, x[1]), ($op).(mlc.edges, x[2]))
        end
    end)
end

"""
    locate(mlc::MultiLeafCollimator, x)

Locate the index `i` such that `mlc.edges[i]<=x<mlc.edges[i+1]`.
"""
locate(mlc::MultiLeafCollimator, x) = locate(mlc.edges, x)

#--- Jaws --------------------------------------------------------------------------------------------------------------

""" 
    Jaws

`Jaws` stores the x and y positions of the jaws.

The x/y positions of the jaws can be accessed through the `getx`/`gety` methods.

The usual constructor directly takes the jaw x and y position vectors, but a
single fieldsize can also be specified.
"""
struct Jaws{T} <: AbstractBeamLimitingDevice
    x::SVector{2, T}
    y::SVector{2, T}

    # Constructors
    function Jaws(x::AbstractVector{T}, y::AbstractVector{T}) where T<:AbstractFloat
        new{T}(SVector{2, T}(x), SVector{2, T}(y))
    end
    function Jaws(x1::T, x2::T, y1::T, y2::T) where T<:AbstractFloat
        new{T}(SVector(x1, x2), SVector(y1, y2))
    end
end

"""
    Jaws(fieldsize::T)

Create jaws with a square field of length `fieldsize`.
"""
Jaws(fieldsize::T) where T<:AbstractFloat = Jaws(-0.5*fieldsize, 0.5*fieldsize, -0.5*fieldsize, 0.5*fieldsize)

getx(jaws::Jaws) = jaws.x
gety(jaws::Jaws) = jaws.y
