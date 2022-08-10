import Base.+, Base.-, Base.getindex, Base.lastindex, Base.size, Base.length
using StaticArrays

export Jaws, MultiLeafCollimator, getx, gety, inaperture, centerposition, edgeposition, extract_subset

"""
    Abstract type for Beam Limiting Devices are based on

`AbstractBeamLimitingDevice` contains the types for MultiLeafCollimators, Jaws, etc.
"""
abstract type AbstractBeamLimitingDevice end

#--- MultiLeafCollimator ------------------------------------------------------

"""
    MultiLeafCollimator

`MultiLeafCollimator` stores the leaf position boundaries (edges) of the MLC.
An MLC consists of a set number of leaf tracks, with two leaf positions for each
track.

The leaf position boundaries are stored in the `positions` vector, which is 1
element longer than the number of leaf tracks, `n`.

Indexing a `MultiLeafCollimator` returns the lower and upper boundaries of the
leaf. Indexing also supports integer ranges. Both support the @views macro.
Iteration goes through each leaf track, returning the lower and upper boundaries.

A `MultiLeafCollimator` can also be shifted through addition or subtraction.
This creates a new `MultiLeafCollimator`.

Other methods include:
- `centerposition`: return the centre of a given leaf track
- `edgeposition`: return the lower edge of a given leaf track
- `width`: return the width of a leaf track
- `locate`: return the leaf index containing the given index
- `extract_subset`: extract a subset of the MLC between an upper and lower bound
"""
struct MultiLeafCollimator{TVector<:AbstractVector} <: AbstractBeamLimitingDevice
    position::TVector
    n::Int
    function MultiLeafCollimator(y)
        n = length(y)-1
        new{typeof(y)}(y, n)
    end
end

# Constructors

"""
    MultiLeafCollimator(n::Int, Δy::Real)

Construct an MLC with `n` number of leaves and of leaf width `Δy`, centered
about zero.
"""
MultiLeafCollimator(n::Int, Δy::Real) = MultiLeafCollimator(Δy*(-0.5*n:0.5*n))

# Size and Length
Base.length(mlc::MultiLeafCollimator) = mlc.n
Base.size(mlc::MultiLeafCollimator) = (length(mlc),)

Base.lastindex(mlc::MultiLeafCollimator) = length(mlc)
Base.eachindex(mlc::MultiLeafCollimator) = Base.OneTo(length(mlc))

Base.show(io::IO, mlc::MultiLeafCollimator) = Base.show(io, mlc.position)

# Indexing
Base.getindex(mlc::MultiLeafCollimator, i::Int) = mlc.position[i:i+1]
Base.getindex(mlc::MultiLeafCollimator, i::UnitRange{Int}) = mlc.position[i[1]:i[end]+1]

# Iteration
Base.iterate(mlc::MultiLeafCollimator) = (mlc[1], 1)
Base.iterate(mlc::MultiLeafCollimator, state) = state>=length(mlc) ? nothing : (mlc[state+1], state+1)

# Views
Base.view(mlc::MultiLeafCollimator, i::UnitRange{Int}) = view(mlc.position, i)
Base.view(mlc::MultiLeafCollimator, i::Int) = view(mlc.position, i:i+1)

# Shifting the MLC
Base.:+(mlc::MultiLeafCollimator, x::Number) = MultiLeafCollimator(mlc.position .+ x)
Base.:-(mlc::MultiLeafCollimator, x::Number) = MultiLeafCollimator(mlc.position .+ x)

# Other Methods
"""
    centerposition(mlc::MultiLeafCollimator, i)

Get the centre of leaf `i`.
"""
centerposition(mlc::MultiLeafCollimator, i::Int) = 0.5*(mlc.position[i] + mlc.position[i+1])

"""
    edgeposition(mlc::MultiLeafCollimator, i)

Get the lower edge of leaf `i`.
"""
edgeposition(mlc::MultiLeafCollimator, i::Int) = mlc.position[i]

"""
    width(mlc::MultiLeafCollimator, i)

Get the width of leaf `i`.
"""
width(mlc::MultiLeafCollimator, i::Int) = mlc.position[i+1] - mlc.position[i]

"""
    locate(mlc::MultiLeafCollimator, x)

Find the index `i` such that `mlc[i][1]<=x<mlc[i][2]`.
"""
locate(mlc::MultiLeafCollimator, x) = locate(mlc.position, x)

"""
    subset_indices(mlc::MultiLeafCollimator, y_lower, y_upper)

Return the indices of the lower/upper leaf track containing `y_lower`/`y_upper`.

These indices can be used to subset the mlc, *e.g* `mlc[i1:i2]`.
Calls `searchsortedlast` on the lower bound `y_lower`, which returns `i` such
that `y[i]<=y<y[i+1]`. Calls `searchsortedfirst` on the upper bound `y_upper`,
which returns `i` such that `y[i-1]<y<=y[i]`.
"""
function subset_indices(mlc::MultiLeafCollimator, y_lower, y_upper)
    i1 = max(1, searchsortedlast(mlc.position, y_lower))
    i2 = min(length(mlc), searchsortedfirst(mlc.position, y_upper)-1)
    i1, i2
end

subset_indices(mlc::MultiLeafCollimator, ylim::AbstractVector) = subset_indices(mlc, ylim...)

"""
    extract_subset(mlc::MultiLeafCollimator, y_lower, y_upper)

Return a subset of the MLC between `y_lower` and `y_upper`.
"""
function extract_subset(mlc::MultiLeafCollimator, y_lower, y_upper)
    i1, i2 = subset_indices(mlc, y_lower, y_upper)
    MultiLeafCollimator(mlc[i1:i2])
end

"""
    extract_subset(args..., ylim::AbstractVector)

Return a subset between `ylim[1]` and `ylim[2]`.
""" extract_subset

extract_subset(mlc::MultiLeafCollimator, ylim::AbstractVector) = extract_subset(mlc, ylim...)

"""
    extract_subset(mlcx, mlc::MultiLeafCollimator, y_lower, y_upper)

Return a subset of the MLC leaf positions (`mlcx`).
"""
function extract_subset(mlcx, mlc::MultiLeafCollimator, y_lower, y_upper)
    i1, i2 = subset_indices(mlc, y_lower, y_upper)
    mlcx[:, i1:i2]
end

extract_subset(mlcx, mlc::MultiLeafCollimator, ylim::AbstractVector) = extract_subset(mlcx, mlc, ylim...)

#--- MultiLeafCollimator Positions -------------------------------------------------------------------------------------

AbstractMLCAperture{T} = AbstractMatrix{T}
AbstractMLCSequence{T} = AbstractArray{T, 3}

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
