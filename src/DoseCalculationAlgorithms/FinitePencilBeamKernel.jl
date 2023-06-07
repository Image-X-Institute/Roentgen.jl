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

"""
    FinitePencilBeamKernel(parameters, scalingfactor, depth, tanθ)

Dose calculation for Finite Pencil Beam Kernel algorithm (Jelen 2005).

Takes the following commissioned parameters:
- `parameters`: Weights and steepness parameters by depth
- `scalingfactor`: Scaling factor matrix by depth and tanθ
- `depth`: Depths of parameters of scaling factor values
- `tanθ`: Angle from central beam axis of scaling factor value
"""
struct FinitePencilBeamKernel{Tparameters<:AbstractInterpolation,
                              TScalingFactor<:AbstractInterpolation,
                              T<:Real} <: AbstractDoseAlgorithm
    parameters::Tparameters
    scalingfactor::TScalingFactor
    α_depth::T  # For scaling depth to interpolation knots
    α_tanθ::T  # For scaling tanθ to interpolation knots
end

function FinitePencilBeamKernel(parameters::AbstractVector{SVector{5, T}},
                                scalingfactor::AbstractMatrix{T},
                                depths::AbstractRange{T}, tanθ::AbstractRange{T}) where T<:Real
    @assert length(depths) == length(parameters)
    @assert length(depths) == size(scalingfactor, 1)
    @assert length(tanθ) == size(scalingfactor, 2)

    interp_method = BSpline(Linear())
    extrap = Interpolations.Line()

    I_paramI = extrapolate(interpolate(parameters, interp_method), extrap)
    I_scalingf = extrapolate(interpolate(scalingfactor, interp_method), extrap)

    α_depth = (length(depths)-1)/depths[end]
    α_tanθ = (length(tanθ)-1)/tanθ[end]

    FinitePencilBeamKernel(I_paramI, I_scalingf, α_depth, α_tanθ)
end

function FinitePencilBeamKernel(parameters::AbstractMatrix, args...)
    @assert size(parameters, 1) == 5
    FinitePencilBeamKernel(SVector{5}.(eachcol(parameters)), args...)
end

Adapt.@adapt_structure FinitePencilBeamKernel

"""
    FinitePencilBeamKernel(filename::String)

Load commissioned parameters from a `.jld` file.
"""
function FinitePencilBeamKernel(filename::String)
    data = JLD2.load(filename)
    FinitePencilBeamKernel(data["parameters"], data["scalingfactor"],
                           data["depth"], data["tantheta"])
end

"""
    calibrate!(calc, MU, fieldsize, SAD[, SSD=SAD])

Calibrate a dose algorithm with given `MU`, `fielsize` and `SAD`.

Scales the dose such that the maximum dose is 1 Gy for `MU` monitor units, given
`fieldsize` and source-axis distance (`SAD`).
Can set source-surface distance `SSD` if `SSD!=SAD`.
"""
function calibrate!(calc::FinitePencilBeamKernel, MU, fieldsize, SAD, SSD=SAD;
                    beamlet_size=5.)

    surf = PlaneSurface(SSD)

    xb = -0.5*fieldsize:beamlet_size:0.5*fieldsize
    bixels = bixel_grid(xb, xb)

    gantry = GantryPosition(0., 0., SAD)

    beamlets = Beamlet.(bixels, (gantry,))

    pos(x) = SVector(0., 0., -x)
    f(x) = -sum(point_dose.(Ref(pos(x)), beamlets, Ref(surf), Ref(calc)))
    result = optimize(f, 0., 100.)
    max_depth_dose = -minimum(result)

    calc.scalingfactor.itp.coefs ./= max_depth_dose*MU
    
    result
end

_scale_clamp(x, α) = max(α*x, 0)+1

function getparams(calc, depth)
    x = _scale_clamp(depth, calc.α_depth)
    p = calc.parameters(x)
    a = p[1]
    ux = SVector(p[2], p[3])
    uy = SVector(p[4], p[5])
    a, ux, uy
end

function getscalingfactor(calc, depth_rad, tanθ)
    x = _scale_clamp(depth_rad, calc.α_depth)
    y = _scale_clamp(abs(tanθ), calc.α_tanθ)
    calc.scalingfactor(x, y)
end

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

    δ = SAD*(r-rₐ)/Rₐ
    x = dot(δ, ax)
    y = dot(δ, ay)

    tanθ = tanθ_to_source(beamlet)

    A = getscalingfactor(calc, depth, tanθ)
    w, ux, uy = getparams(calc, depth)

    F = fpbk_dose(x, y, w, ux, uy, x₀, y₀)
    A*F*(SAD/Rₐ)^2
end
