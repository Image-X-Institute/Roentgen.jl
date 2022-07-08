import Base.+, Base.length, Base.getindex, Base.eachindex
using WriteVTK
using SparseArrays, StaticArrays

export save, DoseGrid, CylinderBounds, gridsize, MeshBounds


# Dose Positions
abstract type DosePositions end

Base.iterate(pos::DosePositions) = iterate(pos, 1)
function Base.iterate(pos::DosePositions, i)
    if(i>length(pos)) return nothing end
    pos[i], i+1
end

function Base.:+(pos::T, x) where T<:DosePositions
    T(pos.x .+ x[1], pos.y .+ x[2], pos.z .+ x[3])
end

# Bounds
struct CylinderBounds{T}
    diameter::T
    height::T
    center::SVector{3, T}
end

function boundingbox(bounds::CylinderBounds)
    xmin, xmax = bounds.center[1] .+ 0.5*bounds.diameter*SVector(-1., 1.)
    ymin, ymax = bounds.center[2] .+ 0.5*bounds.height*SVector(-1., 1.)
    zmin, zmax = bounds.center[3] .+ 0.5*bounds.diameter*SVector(-1., 1.)

    xmin, xmax, ymin, ymax, zmin, zmax
end

function within(bounds::CylinderBounds, p)
    p′ = p - bounds.center
    p′[1]^2 + p′[3]^2 <= 0.25*bounds.diameter^2 && abs(p′[2]) <= 0.5*bounds.height
end

struct MeshBounds{T<:AbstractFloat, TMesh<:SimpleMesh}
    mesh::TMesh
    xlim::SVector{2, T}
    ylim::SVector{2, T}
    zlim::SVector{2, T}
    pad::T
end

function MeshBounds(mesh, pad=10.)
    pos = coordinates.(vertices(mesh))
    x = getindex.(pos, 1)
    y = getindex.(pos, 2)
    z = getindex.(pos, 3)
    xlim = SVector(minimum(x)-pad, maximum(x)+pad)
    ylim = SVector(minimum(y)-pad, maximum(y)+pad)
    zlim = SVector(minimum(z)-pad, maximum(z)+pad)
    MeshBounds(mesh, xlim, ylim, zlim, pad)
end

boundingbox(bounds::MeshBounds) = bounds.xlim..., bounds.ylim..., bounds.zlim...

function within(bounds::MeshBounds, p)
    dmin = minimum(norm.(Ref(Point(p)) .- vertices(bounds.mesh)))
    if(dmin<=bounds.pad)
        return true
    end

    xmin, xmax, ymin, ymax, zmin, zmax = boundingbox(bounds)
    line = Segment(Point(xmin, ymin, zmin), Point(p))
    pI = intersect_mesh(line, bounds.mesh)
    
    length(pI) % 2 != 0
end

# DoseGrid

struct DoseGrid{T}
    x::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
    y::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
    z::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
    indices::Vector{CartesianIndex{3}}
    cells::Vector{Vector{Int}}
end

Base.:+(pos::DoseGrid, x) = DoseGrid(pos.x .+ x[1], pos.y .+ x[2], pos.z .+ x[3], pos.indices)

Base.length(pos::DoseGrid) = length(pos.indices)
Base.size(pos::DoseGrid) = (length(pos),)

gridsize(pos::DoseGrid) = (length(pos.x), length(pos.y), length(pos.z))

Base.getindex(pos::DoseGrid, i::Vararg{Int, 3}) = SVector(pos.x[i[1]], pos.y[i[2]], pos.z[i[3]])
Base.getindex(pos::DoseGrid, i::CartesianIndex{3}) = SVector(pos.x[i[1]], pos.y[i[2]], pos.z[i[3]])
Base.getindex(pos::DoseGrid, i::Int) = pos[pos.indices[i]]

Base.eachindex(pos::DoseGrid) = eachindex(pos.indices)
Base.CartesianIndices(pos::DoseGrid) = pos.indices

function DoseGrid(Δ::T, bounds) where T<:Number
    xmin, xmax, ymin, ymax, zmin, zmax = boundingbox(bounds)

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

    neighbours = vec([CartesianIndex(i,j,k) for i=0:1, j=0:1, k=0:1])[2:end]

    # Cells
    cells = Vector{Int}[]
    for i in eachindex(indices)
        cell = [i]
        for n in neighbours
            index = indices[i] + n

            if(checkbounds(Bool, linear_index, index))
                if(linear_index[index] != 0)
                    push!(cell, linear_index[index])
                end
            end
        end
        if(length(cell)==8)
            push!(cells, cell)
        end
    end

    DoseGrid(x, y, z, indices, cells)
end

function vtk_create_cell(cell)
    n = length(cell)
    c = tuple(cell...)
    if(n==5)
        return MeshCell(VTKCellTypes.VTK_PYRAMID, cell)
    elseif(n==6)
        return MeshCell(VTKCellTypes.VTK_WEDGE, cell)
    elseif(n==8)
        return MeshCell(VTKCellTypes.VTK_VOXEL, cell)
    end
end

function vtk_generate_file(filename, pos::DoseGrid)
    points = [pos[i] for i in eachindex(pos)]
    cells = [vtk_create_cell(cell) for cell in pos.cells]
    vtk_grid(filename, points, cells)
end

function vtk_grid_data(pos::DoseGrid, arr)  
    grid_arr = zeros(size(pos))
    for i in eachindex(arr)
        grid_arr[pos.indices[i]] = arr[i]
    end
    grid_arr
end

function save(filename, pos::DoseGrid, data)
    vtkfile = vtk_generate_file(filename, pos)
    for (key, value) in data
        vtkfile[key, VTKPointData()] = value
    end
    vtk_save(vtkfile)
end
