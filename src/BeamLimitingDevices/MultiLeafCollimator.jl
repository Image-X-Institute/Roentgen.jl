
#--- AbstractMultiLeafCollimator ---------------------------------------------------------------------------------------

abstract type AbstractMultiLeafCollimator <: AbstractBeamLimitingDevice end

export getedges, getpositions, setpositions!, closeleaves!

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

#--- Constructors

"""
    MultiLeafCollimator(n::Int, Δy::Real)

Construct an MLC with `n` number of leaves and of leaf width `Δy`, centered
about zero.
"""
function MultiLeafCollimator(n::Int, Δy::Real)
    x = zeros(2, n)
    y = Δy*(-0.5*n:0.5*n)
    MultiLeafCollimator(x, y)
end

"""
    MultiLeafCollimator(n::Int, Δy::Real)

Construct an MLC with leaf edges.
"""
function MultiLeafCollimator(edges)
    n = length(edges)-1
    x = zeros(2, n)
    MultiLeafCollimator(x, edges)
end

#--- Methods

Base.copy(mlc::MultiLeafCollimator) = MultiLeafCollimator(copy(getpositions(mlc)),
                                                          copy(getedges(mlc)))
Base.similar(mlc::MultiLeafCollimator) = MultiLeafCollimator(getedges(mlc))

# Size and Length
Base.length(mlc::MultiLeafCollimator) = mlc.n
Base.size(mlc::MultiLeafCollimator) = (length(mlc),)

# Iteration
Base.firstindex(mlc::MultiLeafCollimator) = 1
Base.lastindex(mlc::MultiLeafCollimator) = mlc.n
Base.eachindex(mlc::MultiLeafCollimator) = Base.OneTo(mlc.n)

#--- Indexing

Base.getindex(mlc::MultiLeafCollimator, i::Int) = Jaws(mlc.positions[:, i], mlc.edges[i:i+1])

function Base.getindex(mlc::MultiLeafCollimator, i::UnitRange{Int})
    MultiLeafCollimator(mlc.positions[:, i], mlc.edges[i[1]:i[end]+1])
end

function Base.view(mlc::MultiLeafCollimator, i::UnitRange{Int})
    x = @view mlc.positions[:, i]
    y = @view mlc.edges[i[1]:i[end]+1]
    MultiLeafCollimator(x, y)
end

Base.setindex!(mlc::MultiLeafCollimator, x, i) = mlc.positions[:, i] .= x

#--- Methods

locate(mlc::MultiLeafCollimator, x) = locate(mlc.edges, x)

getedges(mlc::MultiLeafCollimator) = mlc.edges
getedges(mlc::MultiLeafCollimator, i::Int) = mlc.edges[i:i+1]
getedges(mlc::MultiLeafCollimator, i::UnitRange{Int}) = mlc.edges[i[1]:i[end]+1]

getpositions(mlc::MultiLeafCollimator) = mlc.positions
getpositions(mlc::MultiLeafCollimator, i) = mlc.positions[:, i]

function setpositions!(mlc::MultiLeafCollimator, x)
    mlc.positions .= x
    mlc
end

closeleaves!(mlc::MultiLeafCollimator) = setpositions!(mlc, zeros(2, length(mlc)))

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

#--- IO

_io_conv(x, lw) = round(Int, lw*(x+120)/240.)
_str_closedaperture(i, lw) = (@sprintf "%3i: " i)*repeat("░", lw-5)

function _str_aperture(i, linepos, lw)
    l = (@sprintf "%3i: " i)*repeat("░", linepos[1]-1)
    c = repeat(" ", linepos[2]-linepos[1]+1)
    r = repeat("░", lw-5-linepos[2])
    l*c*r
end

function Base.show(io::IO, mlc::MultiLeafCollimator)
    println(io, size(mlc.positions, 1), "x", size(mlc.positions, 2), " MultiLeafCollimator")

    maxlines = 20
    linewidth = min(80, displaysize(stdout)[2])

    N = length(mlc)

    truncate_io = N > maxlines

    nlines = min(maxlines, N)

    if truncate_io 
        println(io, "  ⋮")
        indices = max(1, N):min(N, nlines*3÷4)
    else
        indices = eachindex(mlc)
    end

    for i in indices
        pos = mlc.positions[:, i]
        linepos = _io_conv.(pos, linewidth)
        if pos[1] < pos[2]
            println(io, _str_aperture(i, linepos, linewidth))
        else
            println(io, _str_closedaperture(i, linewidth))
        end
    end

    if truncate_io
        println(io, "  ⋮")
    end
end
