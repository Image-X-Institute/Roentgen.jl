module RoentgenCUDAExt

using Roentgen

import Roentgen.AbstractBeamlet, Roentgen.AbstractExternalSurface, Roentgen.AbstractDoseAlgorithm

using CUDA
import Adapt

function dose_fluence_matrix(::Type{CuArray}, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm; kwargs...)
    D = CUDA.zeros(length(pos), length(beamlets))
    dose_fluence_matrix!(D, pos, beamlets, surf, calc; kwargs...)
end

function dose_fluence_matrix!(D::CuArray, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                              maxradius=100.)
    _assert_size(D, pos, beamlets)
    D .= point_dose.(vec(pos), permutedims(vec(beamlets)), Ref(surf), Ref(calc), maxradius)
end

function Adapt.adapt_structure(to, surf::CylindricalSurface)
    cu_rho = Adapt.adapt_structure(to, surf.rho)
    CylindricalSurface(surf.y, surf.ϕ, cu_rho)
end

Adapt.@adapt_structure DoseVolume
Adapt.@adapt_structure LinearSurface
Adapt.@adapt_structure FinitePencilBeamKernel

end
