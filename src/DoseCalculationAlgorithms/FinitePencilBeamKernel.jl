#=  Finite Pencil Beam Kernel

Implements the Jelen et al. 2005, "A finite size pencil beam for IMRT dose optimization" dose calculation kernel.

All references to equations refer to equations in the original paper.
=#

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
