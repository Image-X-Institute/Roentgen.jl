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
