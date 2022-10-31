#=  Finite Pencil Beam Kernel

Implements the Jelen et al. 2005, "A finite size pencil beam for IMRT dose optimization" dose calculation kernel.

All references to equations refer to equations in the original paper.
=#

#--- Abstract Beamlet ---------------------------------------------------------

abstract type AbstractBeamlet <: AbstractFluenceElement end

export Beamlet, FinitePencilBeamKernel

#--- Beamlet ------------------------------------------------------------------

struct Beamlet{T} <: AbstractBeamlet
    halfwidth::SVector{2, T}
    vector::SVector{3, T}
end

"""
    Beamlet(bixel::Bixel, SAD)

Construct a beamlet from a bixel
"""
function Beamlet(bixel::Bixel, SAD)
    halfwidth = 0.5*width(bixel)
    pos = position(bixel)

    a = SVector(pos..., -SAD)
    a = a/norm(a)

    Beamlet(halfwidth, a)
end

struct FinitePencilBeamKernel{TI1, TI2, T} <: AbstractDoseAlgorithm
    w::TI1
    ux::TI2
    uy::TI2
    max_radius::T
end

function FinitePencilBeamKernel(depths, w, ux, uy; max_radius=25.)
    w = LinearInterpolation(depths, w, extrapolation_bc=Interpolations.Line())
    ux = LinearInterpolation(depths, ux, extrapolation_bc=Interpolations.Line())
    uy = LinearInterpolation(depths, uy, extrapolation_bc=Interpolations.Line())
    FinitePencilBeamKernel(w, ux, uy, max_radius)
end

#--- Kernel Profile Functions -------------------------------------------------

@inline fpkb_profile_left(x, u, x₀) = sinh(u*x₀)*exp(u*x)
@inline fpkb_profile_center(x, u, x₀) = 1-exp(-u*x₀)*cosh(u*x)
@inline fpkb_profile_right(x, u, x₀) = sinh(u*x₀)*exp(-u*x)

function fpkb_profile(x, u, x₀)
    x < -x₀ && return fpkb_profile_left(x, u, x₀)
    x > x₀ && return fpkb_profile_right(x, u, x₀)
    return fpkb_profile_center(x, u, x₀)
end

function fpkb_dose(x, y, w, ux, uy, x₀, y₀)
    fx = fpkb_profile.(x, ux, x₀)
    fy = fpkb_profile.(y, uy, y₀)

    w*fx[1]*fy[1] + (1-w)*fx[2]*fy[2]
end

#--- Dose-Fluence Matrix Computations -----------------------------------------

"""
    kernel_size(calc::FinitePencilBeamKernel, pos, beamlets)

Compute the number of bixels with the max. radius of a given position.
"""
function kernel_size(calc::FinitePencilBeamKernel, pos::SVector{3, T}, beamlets, SAD::T) where T<:AbstractFloat
    
    max_kernel_radius² = calc.max_radius^2

    n = 0
    for beamlet in beamlets
        rₐ = dot(pos, beamlet.vector)*beamlet.vector
        R² = dot(pos - rₐ, pos - rₐ)
        if(R² < max_kernel_radius²)
            n += 1
        end
    end
    n
end

"""
    point_kernel!(rowval, nzval, pos::AbstractVector{T}, bixels, surf, calc)

Compute the fluence kernel for a given position.

Designed to be used with a dose-fluence matrix of type SparseMatrixCSC. Stores
the row value in `rowval`, and dose value in `nzval`.
"""
function point_kernel!(rowval, nzval, pos::AbstractVector{T}, beamlets, surf, calc::FinitePencilBeamKernel) where T<:AbstractFloat

    SAD = T(1000)

    max_kernel_radius² = calc.max_radius^2

    depth = getdepth(surf, pos)

    depth < zero(T) && return 0

    w = calc.w(depth)
    ux = calc.ux(depth)
    uy = calc.uy(depth)

    n = 0
    for i in eachindex(beamlets)
        beamlet = beamlets[i]

        Rₐ = dot(pos, beamlet.vector)
        rₐ = Rₐ*beamlet.vector

        x, y = SAD*(pos[1:2] - rₐ[1:2])/Rₐ
        x₀, y₀ = beamlet.halfwidth

        R² = dot(pos - rₐ, pos - rₐ)
        if(R² < max_kernel_radius²)
            n += 1
            rowval[n] = i
            nzval[n] = fpkb_dose(x, y, w, ux, uy, x₀, y₀)/Rₐ^2
        end
    end
    n
end
