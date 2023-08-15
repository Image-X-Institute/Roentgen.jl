abstract type AbstractDoseVolume end

struct DoseVolume{TPos<:AbstractArray, TSurf<:AbstractExternalSurface} <: AbstractDoseVolume
    positions::TPos
    surface::TSurf
end

"""
    getpositions(vol::DoseVolume)

Get dose positions/grid
"""
getpositions(vol::DoseVolume) = vol.positions

"""
    getsurface(vol::DoseVolume)

Get surface
"""
getsurface(vol::DoseVolume) = vol.surface

"""
    dose_fluence_matrix!(D, vol, beamlets, calc)

Gets positions and surface from `vol`
"""
function dose_fluence_matrix!(D, vol::AbstractDoseVolume, beamlets, calc; kwargs...)
    dose_fluence_matrix!(D, getpositions(vol), beamlets, getpositions(surf), calc;
        kwargs...)
end

"""
    dose_fluence_matrix(T, vol, beamlets, calc)

Gets positions and surface from `vol`
"""
function dose_fluence_matrix(T, vol::AbstractDoseVolume, beamlets, calc; kwargs...)
    dose_fluence_matrix(T, getpositions(vol), beamlets, getpositions(surf), calc;
        kwargs...)
end

Adapt.@adapt_structure DoseVolume
