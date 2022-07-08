#
#   Scaled Isoplane Kernel Algorithm
#

using DelimitedFiles, Glob, LsqFit

export ScaledIsoplaneKernel, kernel, norm_depth_dose, calibrate!

mutable struct ScaledIsoplaneKernel{T<:AbstractFloat, TKernel, TPDD} <: AbstractDoseAlgorithm
    kernel::TKernel
    PDD::TPDD
    ref_depth::T
    ref_SSD::T
    max_radius::T
    δsub::SVector{2, T}
end

function ScaledIsoplaneKernel(filename::String, max_kernel_radius; δsub=SVector(1., 1.))

    data = Dict()
    open(filename, "r") do f
        data = JSON.parse(read(f, String))
    end

    ref_depth = data["Ref. Depth"]
    ref_SSD = data["Ref. SSD"]

    # Load Kernel Data
    kernel_radius = data["Kernel"]["Radius"]
    kernel_value = data["Kernel"]["Kernel Value"]

    # Truncate Kernel
    # i = locate(kernel_radius, max_kernel_radius)
    # kernel_radius = kernel_radius[1:i]
    # kernel_value = kernel_value[1:i]
    
    # Calc Log
    @. kernel_value = log(kernel_value)

    # Interpolate onto uniform grid
    Δr = kernel_radius[2]-kernel_radius[1]
    r = kernel_radius[1]:Δr:kernel_radius[end]
    K = LinearInterpolation(kernel_radius, kernel_value)(r)

    # Create interpolation kernel
    kernel = LinearInterpolation(r, K, extrapolation_bc=-Inf)
    
    # Load Depth Dose Data
    pdd_depth = data["Depth Dose"]["Depth"]
    pdd_dose = data["Depth Dose"]["Dose"]
    
    # Log and scale PDD
    @. pdd_dose = log(pdd_dose) - 2*log(100.)

    # Interpolate onto uniform grid
    Δd = pdd_depth[2]-pdd_depth[1]
    d = pdd_depth[1]:Δd:pdd_depth[end]
    D = LinearInterpolation(pdd_depth, pdd_dose)(d)

    # Create interpolation PDD
    PDD = CubicSplineInterpolation(d, D, extrapolation_bc=Interpolations.Line())

    ScaledIsoplaneKernel(kernel, PDD,
                         ref_depth, ref_SSD,
                         max_kernel_radius,
                         δsub)
end

function calibrate!(calc, MU, fieldsize, SSD)
    surf = PlaneSurface(SSD)
    bixel = Bixel(0., fieldsize)
    A = integrate_kernel(calc, bixel, 0., 0.)
    calc.kernel.itp.itp.coefs .-= log(A)
end

"""
    load_kernel_data(filename)

Load a kernel data file in Eclipse format.
"""
function load_kernel_data(filename)
    # Read data file
    data = read_data_file(filename)

    data["radius"], data["value"]
end

"""
    kernel(calc::ScaledIsoplaneKernel, r)

Get the kernel at radius `r`.
"""
kernel(calc::ScaledIsoplaneKernel, r) = exp(calc.kernel(r))

"""
    norm_depth_dose(calc::ScaledIsoplaneKernel, d)

Get the normalised depth dose at depth `d`.
"""
norm_depth_dose(calc::ScaledIsoplaneKernel, d) = exp(calc.PDD(d))

@inline function mayneords_F_factor(SSD, depth, calc)
    ((SSD + calc.ref_depth)*(calc.ref_SSD + depth)/(calc.ref_SSD + calc.ref_depth)/(SSD + depth))^2
end

@inline function inverse_square_correction(SSD, calc)
    ((calc.ref_SSD + calc.ref_depth)/(SSD + calc.ref_depth))^2
end

"""
    kernel_size(calc::ScaledIsoplaneKernel, pos, bixels)

Compute the number of bixels with the max. radius of a given position.
"""
function kernel_size(calc::ScaledIsoplaneKernel, pos::SVector{3, T}, bixels::AbstractVector{<:AbstractBixel{T}}, SAD::T) where T<:AbstractFloat
    
    x_iso, y_iso = scale_to_isoplane(pos, -SAD)

    n = 0
    for i in eachindex(bixels)
        x, y = position(bixels[i])
    
        r² = (x-x_iso)^2 + (y-y_iso)^2
        if(r² < calc.max_radius^2)
            n += 1
        end
    end
    n
end

"""
    subdivide(x, Δx, δxmax)

Subdivide at position `x` with width `Δx` with a max subvision width of `δxmax`

Returns the starting position, the subdivision width and the number of subdivisions.
Allows for use in a loop, where the center of subdivision `i` is at:
    `x0 + i*δx` from `0:nx-1`
"""
@inline function subdivide(x, Δx, δxmax)
    nx = ceil(Int, Δx/δxmax)
    δx = Δx/nx
    x0 = x - (Δx - δx)/2
    x0, δx, nx
end

"""
    integrate_kernel(calc::ScaledIsoplaneKernel, bixel::AbstractBixel{T}, x_iso, y_iso) where T<:AbstractFloat

Integrate the kernel over `bixel` from position `x_iso`, `y_iso`.

Subdivides the bixel for higher accuracy.
"""
function integrate_kernel(calc::ScaledIsoplaneKernel, bixel::AbstractBixel{T}, x_iso::T, y_iso::T) where T<:AbstractFloat
    x0, δx, nx = subdivide(bixel[1], width(bixel, 1), calc.δsub[1])
    y0, δy, ny = subdivide(bixel[2], width(bixel, 2), calc.δsub[2])

    K = zero(T)
    for i=0:nx-1, j=0:ny-1
        r = √((x0 + i*δx - x_iso)^2 + (y0 + j*δy - y_iso)^2)
        K += kernel(calc, r)
    end
    K*δx*δy
end

"""
    point_kernel!(rowval, nzval, pos::AbstractVector{T}, bixels, surf, calc)

Compute the fluence kernel for a given position.

Designed to be used with a dose-fluence matrix of type SparseMatrixCSC. Stores
the row value in `rowval`, and dose value in `nzval`.
"""
function point_kernel!(rowval, nzval, pos::AbstractVector{T}, bixels, surf, calc) where T<:AbstractFloat

    SAD = T(1000.)  # This should be moved, SAD not always == 1000.

    max_kernel_radius² = calc.max_radius^2

    SSD = getSSD(surf, pos)
    depth = getdepth(surf, pos)

    depth < zero(T) && return 0

    PDD = norm_depth_dose(calc, depth)

    cM = mayneords_F_factor(SSD, depth, calc)
    ci = inverse_square_correction(SSD, calc)

    α = PDD*cM*ci
        
    x_iso, y_iso = scale_to_isoplane(pos, -SAD)

    n = 0
    for i in eachindex(bixels)
        bixel = bixels[i]
        x, y = position(bixel)

        r² = (x - x_iso)^2 + (y - y_iso)^2

        if(r² < max_kernel_radius²)
            n += 1
            rowval[n] = i
            nzval[n] = α*integrate_kernel(calc, bixel, x_iso, y_iso)
        end
    end
    n
end
