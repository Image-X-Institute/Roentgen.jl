#=  Finite Pencil Beam Kernel

Implements the Jelen et al. 2005, "A finite size pencil beam for IMRT dose optimization" dose calculation kernel.

All references to equations refer to equations in the original paper.
=#

#--- Abstract beamlet ---------------------------------------------------------

abstract type AbstractBeamlet <: AbstractFluenceElement end

export Beamlet, FinitePencilBeamKernel

#--- Beamlet ------------------------------------------------------------------

struct Beamlet{T} <: AbstractBeamlet
    halfwidth::SVector{2, T}
    ax::SVector{3, T}
    ay::SVector{3, T}
    az::SVector{3, T}
    beamaxis::SVector{3, T}
    SAD::T
    tanθ::T
end

"""
    Beamlet(bixel::Bixel, gantry)

Construct a beamlet from a bixel
"""
function Beamlet(bixel::Bixel, gantry::GantryPosition)

    s = getposition(gantry)
    SAD = getSAD(gantry)

    p = position(bixel)
    hw = 0.5*width(bixel)

    trans = bld_to_fixed(gantry)
    b = trans(SVector(p[1], p[2], -SAD))
    bx = trans(SVector(p[1]+hw[1], p[2], -SAD))
    by = trans(SVector(p[1], p[2]+hw[2], -SAD))

    az = normalize(b - s)
    ax = normalize(cross(az, by-b))
    ay = normalize(cross(az, bx-b))

    tanθ = norm(p)/SAD

    Beamlet(hw, ax, ay, az, s/SAD, SAD, tanθ)
end

source_axis(beamlet::Beamlet) = beamlet.beamaxis
source_axis_distance(beamlet::Beamlet) = beamlet.SAD
beamlet_axes(beamlet::Beamlet) = beamlet.ax, beamlet.ay, beamlet.az

direction(beamlet::Beamlet) = beamlet.az

halfwidth(beamlet::Beamlet) = beamlet.halfwidth
width(beamlet::Beamlet) = 2*halfwidth(beamlet)

source_position(beamlet::Beamlet) = source_axis_distance(beamlet)*source_axis(beamlet)

tanθ_to_source(beamlet::Beamlet) = beamlet.tanθ

struct FinitePencilBeamKernel{Tparameters<:AbstractInterpolation, TScalingFactor<:AbstractInterpolation} <: AbstractDoseAlgorithm
    parameters::Tparameters
    scalingfactor::TScalingFactor
end

function FinitePencilBeamKernel(depths, parameters::AbstractMatrix, tanθ, scalingfactor::AbstractMatrix)
    @assert length(depths) == size(parameters, 2)
    @assert length(depths) == size(scalingfactor, 1)
    @assert length(tanθ) == size(scalingfactor, 2)

    data = SVector{5}.(eachcol(parameters))
    params_interpolator = linear_interpolation(depths, data, extrapolation_bc=Interpolations.Line())

    scalingfactor_interpolator = linear_interpolation((depths, tanθ), scalingfactor, extrapolation_bc=Interpolations.Line())

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

"""
    calibrate!(calc, MU, fieldsize, SAD[, SSD=SAD])

Calibrate a dose algorithm with given `MU`, `fielsize` and `SAD`.

Scales the dose such that the maximum dose is 1 Gy for `MU` monitor units, given
`fieldsize` and source-axis distance (`SAD`).
Can set source-surface distance `SSD` if `SSD!=SAD`.
"""
function calibrate!(calc::FinitePencilBeamKernel, MU, fieldsize, SAD, SSD=SAD; beamlet_size=5.)

    surf = PlaneSurface(SSD)

    xb = -0.5*fieldsize:beamlet_size:0.5*fieldsize
    bixels = bixel_grid(xb, xb)

    gantry = GantryPosition(0., 0., SAD)

    beamlets = Beamlet.(bixels, (gantry,))

    f(x) = -sum(point_dose.(Ref(SVector(0., 0., -x)), beamlets, Ref(surf), Ref(calc)))
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
    ax, ay, a = beamlet_axes(beamlet)

    SAD = source_axis_distance(beamlet)
    s = source_position(beamlet)

    x₀, y₀ = halfwidth(beamlet)

    r = p - s

    Rₐ = dot(r, a)
    rₐ = Rₐ*a
    pₐ = rₐ + s

    depth = getdepth(surf, pₐ, s)
    depth < zero(T) && return zero(T)

    δ = SAD*(r-rₐ)/Rₐ
    x = dot(δ, ax)
    y = dot(δ, ay)

    tanθ = tanθ_to_source(beamlet)

    A = getscalingfactor(calc, depth, tanθ)
    w, ux, uy = getparams(calc, depth)

    F = fpbk_dose(x, y, w, ux, uy, x₀, y₀)
    A*F*(SAD/Rₐ)^2
end
