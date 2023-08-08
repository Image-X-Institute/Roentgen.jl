abstract type AbstractDoseVolume end

getdepth(vol::AbstractDoseVolume) = getdepth.(vol.pos, (vol.surf),)
getSSD(vol::AbstractDoseVolume) = getSSD.(vol.pos, (vol.surf),)

struct DoseVolume{TPos<:AbstractArray, TSurf<:AbstractExternalSurface} <: AbstractDoseVolume
    positions::TPos
    surface::TSurf
end

getpositions(vol::DoseVolume) = vol.positions
getsurface(vol::DoseVolume) = vol.surface

function dose_fluence_matrix!(D, vol::AbstractDoseVolume, beamlets, calc; kwargs...)
    dose_fluence_matrix!(D, getpositions(vol), beamlets, getpositions(surf), calc;
        kwargs...)
end

function dose_fluence_matrix(T, vol::AbstractDoseVolume, beamlets, calc; kwargs...)
    dose_fluence_matrix(T, getpositions(vol), beamlets, getpositions(surf), calc;
        kwargs...)
end

Adapt.@adapt_structure DoseVolume
