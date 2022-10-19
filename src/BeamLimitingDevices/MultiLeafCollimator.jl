
#--- AbstractMultiLeafCollimator ---------------------------------------------------------------------------------------

abstract type AbstractMultiLeafCollimator <: AbstractBeamLimitingDevice end

export getedges, getpositions, setpositions!

# Required for tests.
function Base.:(==)(mlc1::AbstractMultiLeafCollimator, mlc2::AbstractMultiLeafCollimator)
    getpositions(mlc1) == getpositions(mlc2) && getedges(mlc1) == getedges(mlc2)
end

"""
    locate(mlc::AbstractMultiLeafCollimator, x::Number)

Locate the track index `i` which contains the position `x`
""" locate(mlc::AbstractMultiLeafCollimator, x)

#--- MultiLeafCollimator -----------------------------------------------------------------------------------------------

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
struct MultiLeafCollimator{Tpos<:AbstractMatrix, Tedge<:AbstractVector} <: AbstractMultiLeafCollimator
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
MultiLeafCollimator(edges) = MultiLeafCollimator(zeros(2, n), length(edges)-1)

Base.copy(mlc::MultiLeafCollimator) = MultiLeafCollimator(getpositions(mlc), getedges(mlc))
Base.similar(mlc::MultiLeafCollimator) = MultiLeafCollimator(getedges(mlc))

# Size and Length
Base.length(mlc::MultiLeafCollimator) = mlc.n
Base.size(mlc::MultiLeafCollimator) = (length(mlc),)

Base.firstindex(mlc::MultiLeafCollimator) = 1
Base.lastindex(mlc::MultiLeafCollimator) = mlc.n
Base.eachindex(mlc::MultiLeafCollimator) = Base.OneTo(mlc.n)

function show_leaf_positions(io::IO, mlc::MultiLeafCollimator, maxdigits)
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

Base.setindex(mlc::MultiLeafCollimator, x, i::Vararg{2}) = mlc.positions[i] .= x

# Methods

locate(mlc::MultiLeafCollimator, x) = locate(mlc.edges, x)

getedges(mlc::MultiLeafCollimator) = mlc.edges
getedges(mlc::MultiLeafCollimator, i::Int) = mlc.edges[i:i+1]
getedges(mlc::MultiLeafCollimator, i::UnitRange{Int}) = mlc.edges[i[1]:i[end]+1]

getpositions(mlc::MultiLeafCollimator) = mlc.positions
getpositions(mlc::MultiLeafCollimator, i) = mlc.positions[:, i]

setpositions!(mlc::MultiLeafCollimator, x) = vec(mlc.positions) .= vec(x)


function Base.show(io::IO, mlc::MultiLeafCollimator)
    maxdigits = 6
    print(size(mlc.positions, 1), "x", size(mlc.positions, 2), " MultiLeafCollimator")
    show_leaf_positions(io, mlc, maxdigits)
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

