abstract type AbstractDoseVolume end

getdepth(vol::AbstractDoseVolume) = getdepth.(vol.pos, (vol.surf),)
getSSD(vol::AbstractDoseVolume) = getSSD.(vol.pos, (vol.surf),)

struct DoseVolume{TPos<:AbstractArray, TSurf<:AbstractExternalSurface} <: AbstractDoseVolume
    positions::TPos
    surface::TSurf
end

getpositions(vol::DoseVolume) = vol.positions
getsurface(vol::DoseVolume) = vol.surface

