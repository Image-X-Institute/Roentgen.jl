abstract type AbstractDoseVolume end

"""
    DoseVolume

Stores dose positions and an external surface
"""
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
    dose_fluence_matrix!(D, vec(getpositions(vol)), beamlets, getsurface(vol), calc;
        kwargs...)
end

"""
    dose_fluence_matrix(T, vol, beamlets, calc)

Gets positions and surface from `vol`
"""
function dose_fluence_matrix(T, vol::AbstractDoseVolume, beamlets, calc; kwargs...)
    dose_fluence_matrix(T, vec(getpositions(vol)), beamlets, getsurface(vol), calc;
        kwargs...)
end

Adapt.@adapt_structure DoseVolume
