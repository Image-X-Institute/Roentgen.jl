import Base.+, Base.-, Base.getindex, Base.lastindex, Base.size, Base.length
import Base.(==), Base.view

export Jaws, MultiLeafCollimator, getx, gety, inaperture, centerposition, edgeposition, extract_subset

"""
    Abstract type for Beam Limiting Devices are based on

`AbstractBeamLimitingDevice` contains the types for MultiLeafCollimators, Jaws, etc.
"""
abstract type AbstractBeamLimitingDevice end

include("MultiLeafCollimator.jl")
include("MultiLeafCollimatorSequence.jl")

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
