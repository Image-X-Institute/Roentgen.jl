import Base.+, Base.length, Base.getindex, Base.eachindex
using WriteVTK
using SparseArrays, StaticArrays

export save, DoseGrid, DoseGridMasked, CylinderBounds, gridsize, MeshBounds

#--- Bounds ------------------------------------------------------------------------------------------------------------

abstract type AbstractBounds end

"""
    within(bounds::AbstractBounds, p)

Returns `true` if the point `p` is within bounds
""" within(bounds::AbstractBounds, p)

"""
    boundsextent(bounds::AbstractBounds)

Returns the bounding box of the bounds as a 6 element tuple:
`xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`
""" boundsextent(bounds::AbstractBounds)

"""
    CylinderBounds{T}

Represents a cylinder about the `z` axis.
"""
struct CylinderBounds{T} <: AbstractBounds
    diameter::T
    height::T
    center::SVector{3, T}

    function CylinderBounds(diameter, height)
        T = typeof(diameter)
        new{T}(diameter, height, zeros(SVector{3}))
    end

    function CylinderBounds(diameter, height, center)
        T = typeof(diameter)
        new{T}(diameter, height, SVector{3}(center))
    end
end

"""
    boundsextent(bounds::CylinderBounds)

For a cylinder
"""
function boundsextent(bounds::CylinderBounds)
    xmin, xmax = bounds.center[1] .+ 0.5*bounds.diameter*SVector(-1., 1.)
    ymin, ymax = bounds.center[2] .+ 0.5*bounds.diameter*SVector(-1., 1.)
    zmin, zmax = bounds.center[3] .+ 0.5*bounds.height*SVector(-1., 1.)

    xmin, xmax, ymin, ymax, zmin, zmax
end

"""
    within(bounds::CylinderBounds, p)

Whether `p` is within the cylinder
"""
function within(bounds::CylinderBounds, p)
    p′ = p - bounds.center
    p′[1]^2 + p′[2]^2 <= 0.25*bounds.diameter^2 && abs(p′[3]) <= 0.5*bounds.height
end

"""
    MeshBounds{T, TMesh}

Represents a mesh surface.
"""
struct MeshBounds{T<:AbstractFloat, TMesh<:SimpleMesh} <: AbstractBounds
    mesh::TMesh
    box::SVector{6, T}
    pad::T
end

"""
    MeshBounds

Construct a `MeshBounds` with a mesh

Optionally, can specify a `pad` which defines the distance from the mesh that is
still considered within bounds.
"""
function MeshBounds(mesh, pad=10.) <: AbstractBounds
    pos = coordinates.(vertices(mesh))
    x = getindex.(pos, 1)
    y = getindex.(pos, 2)
    z = getindex.(pos, 3)
    box = SVector(minimum(x)-pad, maximum(x)+pad,
                  minimum(y)-pad, maximum(y)+pad,
                  minimum(z)-pad, maximum(z)+pad)
    MeshBounds(mesh, box, pad)
end

"""
    boundsextent(bounds::MeshBounds)

For a mesh
"""
boundsextent(bounds::MeshBounds) = bounds.box

"""
    within(bounds::MeshBounds, p)

Whether `p` is within the mesh
"""
function within(bounds::MeshBounds, p)
    dmin = minimum(norm.(Ref(Point(p)) .- vertices(bounds.mesh)))
    dmin<=bounds.pad && return true

    xmin, xmax, ymin, ymax, zmin, zmax = boundsextent(bounds)
    line = Segment(Point(xmin, ymin, zmin), Point(p))
    pI = intersect_mesh(line, bounds.mesh)
    
    length(pI) % 2 != 0
end

#--- Dose Positions ----------------------------------------------------------------------------------------------------

abstract type DosePositions end

Base.iterate(pos::DosePositions) = iterate(pos, 1)
function Base.iterate(pos::DosePositions, i)
    i>length(pos) && return nothing
    pos[i], i+1
end

#--- AbstractDoseGrid ----------------------------------------------------------------------------------------------------------

"""
    AbstractDoseGrid <: DosePositions

A type of Dose Positions on a regular grid.

Must have the fields `x`, `y`, and `z`.
"""
abstract type AbstractDoseGrid <: DosePositions end

Base.getindex(pos::AbstractDoseGrid, i::Vararg{Int, 3}) = SVector(getindex.(pos.axes, i)...)
Base.getindex(pos::AbstractDoseGrid, i::CartesianIndex{3}) = pos[i[1], i[2], i[3]]

#--- DoseGrid ----------------------------------------------------------------------------------------------------------

struct DoseGrid{TVec<:AbstractVector} <: AbstractDoseGrid
    axes::SVector{3, TVec}
    DoseGrid(x, y, z) = new{typeof(x)}(SVector(x, y, z))
end

function DoseGrid(Δ, bounds::AbstractBounds)
    xmin, xmax, ymin, ymax, zmin, zmax = boundsextent(bounds)

    x = xmin:Δ:xmax
    y = ymin:Δ:ymax
    z = zmin:Δ:zmax

    DoseGrid(x, y, z)
end

Base.size(pos::DoseGrid) = tuple(length.(pos.axes)...)
Base.length(pos::DoseGrid) = prod(length.(pos.axes))

Base.eachindex(pos::DoseGrid) = Base.OneTo(length(pos))
Base.CartesianIndices(pos::DoseGrid) = CartesianIndices(size(pos))

Base.getindex(pos::DoseGrid, i::Int) = pos[CartesianIndices(pos)[i]]

#-- IO 

function save(filename::String, pos::DoseGrid, data::Vararg)
    vtk_grid(filename, pos.axes...) do vtkfile
        for (key, value) in data
            vtkfile[key, VTKPointData()] = value
        end
    end
end

#--- DoseGridMasked ----------------------------------------------------------------------------------------------------

struct DoseGridMasked{TVec<:AbstractVector} <: AbstractDoseGrid
    axes::SVector{3, TVec}
    indices::Vector{CartesianIndex{3}}
    cells::Vector{Vector{Int}}
end

function DoseGridMasked(Δ, bounds::AbstractBounds)
    xmin, xmax, ymin, ymax, zmin, zmax = boundsextent(bounds)

    x = xmin:Δ:xmax
    y = ymin:Δ:ymax
    z = zmin:Δ:zmax

    grid_index = CartesianIndices( (length(x), length(y), length(z)) )
    indices = CartesianIndex{3}[]
    linear_index = zeros(Int, length(x), length(y), length(z))

    # Points
    for i in grid_index
        p = SVector(x[i[1]], y[i[2]], z[i[3]])
        if(within(bounds, p))
            push!(indices, grid_index[i])
            linear_index[i] = length(indices)
        end
    end

    # Cells

    neighbours = vec([CartesianIndex(i,j,k) for i=0:1, j=0:1, k=0:1])[2:end]

    cells = Vector{Int}[]
    for i in eachindex(indices)
        cell = [i]
        for n in neighbours
            index = indices[i] + n

            if(checkbounds(Bool, linear_index, index) && linear_index[index] != 0)
                push!(cell, linear_index[index])
            end
        end
        if(length(cell)==8)
            push!(cells, cell)
        end
    end

    DoseGridMasked(SVector(x, y, z), indices, cells)
end

Base.:+(pos::DoseGridMasked, x) = DoseGridMasked(pos.x .+ x[1], pos.y .+ x[2], pos.z .+ x[3], pos.indices)

Base.length(pos::DoseGridMasked) = length(pos.indices)
Base.size(pos::DoseGridMasked) = (length(pos),)

gridsize(pos::DoseGridMasked) = (length(pos.x), length(pos.y), length(pos.z))

Base.getindex(pos::DoseGridMasked, i::Int) = pos[pos.indices[i]]

Base.eachindex(pos::DoseGridMasked) = eachindex(pos.indices)
Base.CartesianIndices(pos::DoseGridMasked) = pos.indices

#--- IO

function Base.show(io::IO, pos::DoseGridMasked)
    print(io, "x=")
    Base.show(io, pos.axes[1])
    print(io, " y=")
    Base.show(io, pos.axes[2])
    print(io, " z=")
    Base.show(io, pos.axes[3])
    println(io, " npts=", length(pos.indices))
end

function vtk_create_cell(cell)
    n = length(cell)
    n==5 && return MeshCell(VTKCellTypes.VTK_PYRAMID, cell)
    n==6 && return MeshCell(VTKCellTypes.VTK_WEDGE, cell)
    n==8 && return MeshCell(VTKCellTypes.VTK_VOXEL, cell)
end

function save(filename::String, pos::DoseGridMasked, data::Vararg)
    points = [pos[i] for i in eachindex(pos)]
    cells = vtk_create_cell.(pos.cells)

    vtk_grid(filename, points, cells) do vtkfile
        for (key, value) in data
            vtkfile[key, VTKPointData()] = value
        end
    end
end
save(filename::String, pos::DoseGridMasked, data::Dict) =  save(filename, pos, data...)
