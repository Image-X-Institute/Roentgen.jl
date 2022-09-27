export MultiLeafCollimatorSequence

"""
    MultiLeafCollimatorSequence

`MultiLeafCollimatorSequence` stores the leaf positions and edges of an MLC.

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
struct MultiLeafCollimatorSequence{Tpos<:AbstractArray, Tedge<:AbstractVector} <: AbstractMultiLeafCollimator
    positions::Tpos
    edges::Tedge
    nleaves::Int
    napertures::Int
    function MultiLeafCollimatorSequence(positions::AbstractArray, edges::AbstractVector)
        _, nleaves, napertures = size(positions)
        @assert size(positions, 1) == 2 "Number of leaf positions per track != 2"
        @assert size(positions, 2) == length(edges)-1 "Length of positions and edges do not match"
        new{typeof(positions), typeof(edges)}(positions, edges, nleaves, napertures)
    end
end

function MultiLeafCollimatorSequence(napertures::Int, edges::AbstractVector)
    nleaves = length(edges)-1
    positions = zeros(2, nleaves, napertures)
    MultiLeafCollimatorSequence(positions, edges)
end

# Size and Length
Base.length(mlc::MultiLeafCollimatorSequence) = mlc.napertures
Base.size(mlc::MultiLeafCollimatorSequence) = (length(mlc),)

Base.firstindex(mlc::MultiLeafCollimatorSequence) = 1
Base.lastindex(mlc::MultiLeafCollimatorSequence) = mlc.napertures
Base.eachindex(mlc::MultiLeafCollimatorSequence) = Base.OneTo(mlc.napertures)

# Indexing
Base.getindex(mlc::MultiLeafCollimatorSequence, i::Int) = MultiLeafCollimator(mlc.positions[:, :, i], mlc.edges)
Base.getindex(mlc::MultiLeafCollimatorSequence, i::UnitRange{Int}) = MultiLeafCollimatorSequence(mlc.positions[:, :, i], mlc.edges)

# Views
Base.view(mlc::MultiLeafCollimatorSequence, i::Int) = MultiLeafCollimator((@view mlc.positions[:, :, i]), mlc.edges)
Base.view(mlc::MultiLeafCollimatorSequence, i::UnitRange{Int}) = MultiLeafCollimatorSequence((@view mlc.positions[:, :, i]), mlc.edges)

function Base.show(io::IO, mlc::MultiLeafCollimatorSequence)
    maxdigits = 6
    println(join(size(mlc.positions), "x", "x"), " MultiLeafCollimatorSequence")
    if(length(mlc)<=5)
        for n in eachindex(mlc)
            print(io, "n = $n")
            show_leaf_positions(io, mlc[n], maxdigits)
            println(io)
        end
    else
        for n in eachindex(mlc[1:2])
            print(io, "n = $n")
            show_leaf_positions(io, mlc[n], maxdigits)
            println(io)
        end
        for _=1:3
            println(io, "⋮")
        end
        for n in eachindex(mlc[end-1:end])
            print(io, "n = $(lastindex(mlc)-2+n)")
            show_leaf_positions(io, mlc[n], maxdigits)
            println(io)
        end
    end
end

# Operations
for op in (:+, :-, :*, :/)
    eval(quote
        function Base.$op(mlc::MultiLeafCollimatorSequence, x::Union{AbstractVector,Tuple})
            MultiLeafCollimatorSequence( ($op).(mlc.positions, x[1]), ($op).(mlc.edges, x[2]))
        end
    end)
end
