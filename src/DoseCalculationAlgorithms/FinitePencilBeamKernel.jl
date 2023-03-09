#=  Finite Pencil Beam Kernel

Implements the Jelen et al. 2005, "A finite size pencil beam for IMRT dose optimization" dose calculation kernel.

All references to equations refer to equations in the original paper.
=#

#--- Abstract beamlet ---------------------------------------------------------

abstract type Abstractbeamlet <: AbstractFluenceElement end

export Beamlet, FinitePencilBeamKernel

#--- Beamlet ------------------------------------------------------------------

struct Beamlet{T} <: Abstractbeamlet
    halfwidth::SVector{2, T}
    vector::SVector{3, T}
    position::SVector{3, T}
end

"""
    Beamlet(bixel::Bixel, gantry)

Construct a beamlet from a bixel
"""
function Beamlet(bixel::Bixel, gantry::GantryPosition)
    SAD = getSAD(gantry)

    halfwidth = 0.5*width(bixel)

    s = getposition(gantry)
    b = bld_to_fixed(gantry)(SVector(position(bixel)..., -SAD))

    a = normalize(b - s)

    Beamlet(halfwidth, a, s)
end

getposition(beamlet::Beamlet) = beamlet.position
getdirection(beamlet::Beamlet) = beamlet.vector
gethalfwidth(beamlet::Beamlet) = beamlet.halfwidth
getwidth(beamlet::Beamlet) = 2*gethalfwidth(beamlet)

struct FinitePencilBeamKernel{Tparameters<:AbstractInterpolation, TScalingFactor<:AbstractInterpolation} <: AbstractDoseAlgorithm
    parameters::Tparameters
    scalingfactor::TScalingFactor
end

function FinitePencilBeamKernel(depths, parameters::AbstractMatrix, tanθ, scalingfactor::AbstractMatrix)
    @assert length(depths) == size(parameters, 2)
    @assert length(depths) == size(scalingfactor, 1)
    @assert length(tanθ) == size(scalingfactor, 2)

    data = SVector{5}.(eachcol(parameters))
    params_interpolator = LinearInterpolation(depths, data, extrapolation_bc=Interpolations.Line())

    scalingfactor_interpolator = LinearInterpolation((depths, tanθ), scalingfactor, extrapolation_bc=Interpolations.Line())

    FinitePencilBeamKernel(params_interpolator, scalingfactor_interpolator)
end

function FinitePencilBeamKernel(fid::HDF5.H5DataStore)
    parameters = read(fid["parameters"])
    depths = read(fid["depth"])

    scalingfactor = read(fid["scaling_factor"])
    tanθ = read(fid["tan_theta"])

    FinitePencilBeamKernel(depths, parameters, tanθ, scalingfactor)
end

function FinitePencilBeamKernel(filename::String; fieldsize=100.)
    h5open(filename, "r") do fid
        dset = fid["fieldsize-$(Int(fieldsize))mm"]
        FinitePencilBeamKernel(dset)
    end
end

function calibrate!(calc::FinitePencilBeamKernel, MU, fieldsize, SAD, SSD=SAD)

    surf = PlaneSurface(SSD)

    hw = 0.5*fieldsize*SVector(1., 1.)
    a = SVector(0., 0., -1.)
    s = SVector(0., 0., SAD)
    beamlet = Beamlet(hw, a, s)

    f(x) = -point_dose(SVector(0., 0., -x), beamlet, surf, calc)
    result = optimize(f, 0., 100.)
    max_depth_dose = -minimum(result)

    calc.scalingfactor.itp.coefs ./= max_depth_dose*MU
    
    result
end

function getparams(calc, depth)
    p = calc.parameters(depth)
    a = p[1]
    ux = SVector(p[2], p[3])
    uy = SVector(p[4], p[5])
    a, ux, uy
end

getscalingfactor(calc, depth_rad, tanθ) = calc.scalingfactor(depth_rad, abs(tanθ))

#--- Kernel Profile Functions -------------------------------------------------

@inline fpbk_profile_left(x, u, x₀) = sinh(u*x₀)*exp(u*x)
@inline fpbk_profile_center(x, u, x₀) = 1-exp(-u*x₀)*cosh(u*x)
@inline fpbk_profile_right(x, u, x₀) = sinh(u*x₀)*exp(-u*x)

function fpbk_profile(x, u, x₀)
    x < -x₀ && return fpbk_profile_left(x, u, x₀)
    x > x₀ && return fpbk_profile_right(x, u, x₀)
    return fpbk_profile_center(x, u, x₀)
end

function fpbk_dose(x, y, w, ux, uy, x₀, y₀)
    fx = fpbk_profile.(x, ux, x₀)
    fy = fpbk_profile.(y, uy, y₀)

    w*fx[1]*fy[1] + (1-w)*fx[2]*fy[2]
end


function point_dose(p::SVector{3, T}, beamlet::Beamlet, surf::AbstractExternalSurface, calc::FinitePencilBeamKernel) where T<:AbstractFloat
    a = getdirection(beamlet)
    s = getposition(beamlet)
    x₀, y₀ = gethalfwidth(beamlet)

    SAD = norm(s)
    r = p - s

    Rₐ = dot(r, a)
    rₐ = Rₐ*a + s

    depth = getdepth(surf, rₐ, s)
    depth < zero(T) && return zero(T)

    sy = SVector(0, 1, 0)
    sx = cross(sy, s/SAD)

    δ = SAD*(r - Rₐ*a)/Rₐ
    x, y = dot(δ, sx), dot(δ, sy)

    cosθ = dot(a, s/SAD)
    tanθ = √(1-cosθ^2)/cosθ

    w, ux, uy = getparams(calc, depth)
    A = getscalingfactor(calc, depth, tanθ)
    
    A*fpbk_dose(x, y, w, ux, uy, x₀, y₀)*(SAD/Rₐ)^2
end
