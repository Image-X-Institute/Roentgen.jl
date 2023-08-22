"""
    MockKernel

A mock dose algorithm used as an example and tests
"""
struct MockKernel <: AbstractDoseAlgorithm end

calibrate!(calc::MockKernel, args...) = nothing

function point_dose(p::SVector{3, T}, beamlet, surf, calc::MockKernel) where T<:Real
    SAD = source_axis_distance(beamlet)
    s = source_position(beamlet)
    ax, ay, a = beamlet_axes(beamlet)

    x₀, y₀ = halfwidth(beamlet)

    r = p - s

    Rₐ = dot(r, a)
    rₐ = Rₐ*a

    δ = SAD*(r-rₐ)/Rₐ
    x = dot(δ, ax)
    y = dot(δ, ay)

    abs(x)<x₀ && abs(y)<y₀ && return one(T)
    T(1e-3)
end
