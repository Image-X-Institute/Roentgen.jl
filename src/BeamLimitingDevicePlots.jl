""" Plot Functions for Beam Limiting Devices

Plots jaw positions and multileafcollimator apertures.
Sets axes limits to jaw positions
"""

export plot_bld, plot_bld!, axes_lims!, axes_lims

using Plots

#--- Generic BeamLimitingDevice Plot -----------------------------------------------------------------------------------
"""
    plot_bld

Main plot command for beam limiting devices (e.g. jaws, MLCs).

See ?plot_bld for more information.
"""
plot_bld(args...; kwargs...) = plot_bld!(plot(), args...; kwargs...)

"""
    plot_bld!

Main plot command for beam limiting devices (e.g. jaws, MLCs).

Creates a 2D plot of the given beam limiting device in the beam limiting device
(BLD) coordinate system, showing their relevant positions (e.g. aperture
positions for MLC)
    
Use `plot_bld` to create a new plot object, and `plot_bld!` to add to an
existing one.
"""
plot_bld!(args...; kwargs...) = plot_bld!(plot!(), args...; kwargs...)


#--- Jaws Plot ---------------------------------------------------------------------------------------------------------

"""
    plot_bld!(p, jaws::Jaws; kwargs...)

When applied to jaws, plots a box indicating the jaw positions
"""
function plot_bld!(p, jaws::Jaws; kwargs...)
    # Create a box out of the four jaw positions
    x = vcat(jaws.x[1], jaws.x[2], jaws.x[2], jaws.x[1], jaws.x[1])
    y = vcat(jaws.y[1], jaws.y[1], jaws.y[2], jaws.y[2], jaws.y[1])

    plot!(p, x, y; aspect_ratio=1, kwargs...)
end

#--- MultiLeafCollimator Plot ------------------------------------------------------------------------------------------

"""
    plot_bld!(p, mlcx, mlc::MultiLeafCollimator; kwargs...)

When applied to MultiLeafCollimator with leaf positions, plot the open aperture.

Fills the area obstructed by leaves, unless `invert=true` where it fills the open aperture.
"""
function plot_bld!(p, mlcx, mlc::MultiLeafCollimator; invert=false, leaf_length=125., fill=true, fillalpha=0.1, kwargs...)

    # Set fill=false if fillalpha = 0. (completely transparent)
    if(fill)
        if(fillalpha>0.)
            fill=true
        else
            fill=false
        end
    end

    # Segment the MLC positions
    x, y = segment_mlc(mlcx, mlc)

    # Plot the first segment
    plot!(p, x[1], y[1]; aspect_ratio=1, kwargs...)
     # If there are more segments, plot them too. primary=False ensures same style as first
    if(length(x)>1)
        plot!(p, x[2:end], y[2:end]; primary=false, kwargs...)
    end

    # If invert==true, fill the inside of the aperture
    if(invert)
        plot!(p, x, y; line=false, primary=false, fill=fill, fillalpha=fillalpha, kwargs...)
    # Otherwise, fill the outside, up to the leaf length
    else
        # Plot the outside edge of the mlc
        x, y = segment_mlc(mlcx .+ leaf_length*[-1., 1.], mlc)
        plot!(p, x, y; primary=false, kwargs...)
        
        # Fill the areas obstructed by the leaves
        if(fill)
            xB = segmentx(mlcx[1, :], mlcx[1, :].-leaf_length)
            xA = segmentx(mlcx[2, :], mlcx[2, :].+leaf_length)
            y = segmenty(mlc[1:end])
            plot!(p, xB, y; line=false, primary=false, fill=true, fillalpha=fillalpha, kwargs...)
            plot!(p, xA, y; line=false, primary=false, fill=true, fillalpha=fillalpha, kwargs...)
        end
    end
    p
end

"""
    segmentx(xL, xR)

Creates a vector of MLC x positions, used for plotting the position of the
aperture.
"""
function segmentx(xL::Vector{T}, xR::Vector{T}) where T<:Number
    xs = T[]
    n = length(xL)
    for i=1:n
        push!(xs, xL[i])
        push!(xs, xL[i])
    end
    for i=n:-1:1
        push!(xs, xR[i])
        push!(xs, xR[i])
    end
    push!(xs, xL[1])
    xs
end

"""
    segmenty(y)

Creates a vector of MLC y positions, used for plotting the position of the
aperture.
"""
function segmenty(y::AbstractVector{T}) where T<:Number
    ys = T[]
    n = length(y)
    for i=1:n-1
        push!(ys, y[i])
        push!(ys, y[i+1])
    end
    for i=n-1:-1:1
        push!(ys, y[i+1])
        push!(ys, y[i])
    end
    push!(ys, y[1])
    ys
end

"""
    segment_indices(mlcx)

Segment the MLC aperture into separate contained "openings", returning the leaf
index.

The open aperture separates when xB[i] > xA[i+1] or xA[i] < xB[i+1]. This index
is used to remove plotting artefacts.

# Examples
In this case, the aperture has two "openings", at x=-1.->1. and at x=2.->4.
These openings are separated in the y direction, but this information is not
required to separate the openings. 
```julia-repl
julia> mlcx = [-1. -1. 2. 2.;
                1.  1. 4. 4.];
julia> segment_indices(mlcx)
3-element Vector{Int64}:
 1
 3
 5
```
The first aperture opening runs from index 1 -> 2, and the second opening
from index 3 -> 4
"""
function segment_indices(mlcx)
    indices = [1]
    n = size(mlcx, 2)
    for i=1:n-1
        if(mlcx[1, i+1] > mlcx[2, i] || mlcx[1, i] > mlcx[2, i+1])
            push!(indices, i+1)
        end
    end
    push!(indices, n+1)
    indices
end

"""
    segment_mlc(mlcx, mlc)

Segment the MLC aperture into separate contained "openings", return the x and y
positions of the aperture.

Calls `segment_indices` to get the indices for each section, then `segmentx` 
and `segmenty` to retrieve the positions of the aperture.
"""
function segment_mlc(mlcx::AbstractMatrix{T}, mlc) where T<:AbstractFloat
    indices = segment_indices(mlcx)
    
    x = Vector{T}[]
    y = Vector{T}[]

    for i=1:length(indices)-1
        i_start = indices[i]
        i_end = indices[i+1]-1

        xi = segmentx(mlcx[1, i_start:i_end], mlcx[2, i_start:i_end])
        yi = segmenty(mlc[i_start:i_end])
        
        push!(x, xi)
        push!(y, yi)
    end
    x, y
end

"""
    plot_bld!(p, mlc::MultiLeafCollimator; kwargs...)

When applied to MultiLeafCollimator, plot a series of horizontal lines at each leaf y edge.
"""
plot_bld!(p::Plots.Plot, mlc::MultiLeafCollimator; kwargs...) = hline!(p, gety(mlc); kwargs...)


#--- Axes Limits -------------------------------------------------------------------------------------------------------

"""
    axes_lims!(p, jaws::Jaws; pad=10.)

Set the axes limits to the position of the jaws.

By default, a pad of 10 mm is added to each side.
"""
function axes_lims!(p, jaws::Jaws; pad=10.)
    padding = pad.*[-1, 1]
    plot!(p, xlim=jaws.x .+ padding, ylim=jaws.y .+ padding)
end

"""
    axes_lims!(jaws::Jaws; pad=10.)

When no plot object given, applies to latest plot created.
"""
axes_lims!(jaws::Jaws; pad=10.) = axes_lims!(plot!(), jaws; pad=pad)

#--- Fluence Grid ------------------------------------------------------------------------------------------------------

function edgepositions(bixels, i)
    x = position.(bixels, i)
    w = width.(bixels, i)
    vcat(x[1]-0.5*w[1], x .+ 0.5*w)
end

function plot_bld!(p, bixelgrid::AbstractMatrix{<:AbstractBixel}, Ψ::AbstractMatrix; kwargs...)
    x = edgepositions((@view bixelgrid[:, 1]), 1)
    y = edgepositions((@view bixelgrid[1, :]), 2)

    heatmap!(p, x, y, Ψ'; kwargs...)
end