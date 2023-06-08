import Base.+, Base.length, Base.getindex, Base.eachindex

export save, DoseGrid, DoseGridMasked, CylinderBounds, gridsize, MeshBounds, SurfaceBounds

export getaxes

#--- Bounds ------------------------------------------------------------------------------------------------------------

abstract type AbstractBounds end

"""
    within(bounds::AbstractBounds, p)

Returns `true` if the point `p` is within bounds
""" within(bounds::AbstractBounds, p)

"""
    extent(bounds::AbstractBounds)

Returns the bounding box of the bounds as two vectors:
`pmin`, `pmax`
""" extent(bounds::AbstractBounds)

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
    extent(bounds::CylinderBounds)

For a cylinder
"""
function extent(bounds::CylinderBounds)
    xmin, xmax = bounds.center[1] .+ 0.5*bounds.diameter*SVector(-1., 1.)
    ymin, ymax = bounds.center[2] .+ 0.5*bounds.diameter*SVector(-1., 1.)
    zmin, zmax = bounds.center[3] .+ 0.5*bounds.height*SVector(-1., 1.)

    SVector(xmin, ymin, zmin), SVector(xmax, ymax, zmax)
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

Represents a surface from a mesh.
"""
struct MeshBounds{T<:AbstractFloat, TMesh<:SimpleMesh} <: AbstractBounds
    mesh::TMesh
    box::SVector{6, T}
    pad::T

    function MeshBounds(mesh, pad=0.)
        pos = coordinates.(vertices(mesh))
        x = getindex.(pos, 1)
        y = getindex.(pos, 2)
        z = getindex.(pos, 3)
        box = SVector(minimum(x)-pad, maximum(x)+pad,
                      minimum(y)-pad, maximum(y)+pad,
                      minimum(z)-pad, maximum(z)+pad)

        T = eltype(pos[1])
        Tmesh = typeof(mesh)
        new{T, Tmesh}(mesh, box, pad)
    end

end

"""
    extent(bounds::MeshBounds)

For a mesh
"""
function extent(bounds::MeshBounds)
    # box = boundingbox(bounds.mesh)
    # coordinates(minimum(box)), coordinates(maximum(box))
    idx = SVector(1, 3, 5)
    bounds.box[idx], bounds.box[idx.+1]
end

"""
    within(bounds::MeshBounds, p)

Whether `p` is within the mesh
"""
function within(bounds::MeshBounds, p)
    dmin = minimum(norm.(Ref(Point(p)) .- vertices(bounds.mesh)))
    dmin<=bounds.pad && return true

    pmin, _ = extent(bounds)
    line = Segment(Point(pmin), Point(p))
    pI = intersect_mesh(line, bounds.mesh)
    
    length(pI) % 2 != 0
end


"""
    SurfaceBounds{TSurface}

Represents a surface from a mesh.
"""
struct SurfaceBounds{TSurface<:AbstractExternalSurface} <: AbstractBounds
    surf::TSurface
end

"""
    extent(bounds::SurfaceBounds)

For a mesh
"""
extent(bounds::SurfaceBounds) = extent(bounds.surf)

"""
    within(bounds::SurfaceBounds, p)

Whether `p` is within the mesh
"""
within(bounds::SurfaceBounds, p) = isinside(bounds.surf, p)

#--- Dose Positions ----------------------------------------------------------------------------------------------------

abstract type DosePositions end

Base.iterate(pos::DosePositions) = iterate(pos, 1)
function Base.iterate(pos::DosePositions, i)
    i>length(pos) && return nothing
    pos[i], i+1
end

for op in (:+, :-)
    eval(quote
        function Base.$op(pos::DosePositions, v::AbstractVector{<:AbstractVector})
            ($op).(pos, v)
        end
        function Base.$op(pos::DosePositions, v::AbstractVector{<:Real})
            ($op).(pos, (v,))
        end
    end)
end

#--- AbstractDoseGrid ----------------------------------------------------------------------------------------------------------

"""
    AbstractDoseGrid <: DosePositions

A type of Dose Positions on a regular grid.

Must have the fields `x`, `y`, and `z`.
"""
abstract type AbstractDoseGrid <: DosePositions end

Base.getindex(pos::AbstractDoseGrid, i::Vararg{Int, 3}) = SVector(getindex.(pos.axes, i))
Base.getindex(pos::AbstractDoseGrid, i::CartesianIndex{3}) = pos[i[1], i[2], i[3]]

"""
    getaxes(pos::AbstractDoseGrid[, dim])

Return the axes of the grid.

Optionally, can specify which dimension.
"""
getaxes(pos::AbstractDoseGrid) = pos.axes
getaxes(pos::AbstractDoseGrid, dim) = pos.axes[dim]

#--- DoseGrid ----------------------------------------------------------------------------------------------------------
"""
    DoseGrid

Cartesian Dose Grid
"""
struct DoseGrid{TVec<:AbstractVector} <: AbstractDoseGrid
    axes::SVector{3, TVec}
    function DoseGrid(x, y, z)
        ax = x, y, z
        if(!all(@. typeof(ax) <: StepRangeLen))
            ax = collect.(ax)
        end
        new{typeof(ax[1])}(SVector(ax...))
    end
end

"""
    DoseGrid(Δ, bounds::AbstractBounds)

Construct a `DoseGrid` with spacing `Δ` within `bounds`.

Can specify a coordinate system transformation in `transform` to generate dose
points in a different coordinate system to the bounds.
"""
function DoseGrid(Δ, bounds::AbstractBounds, transform=IdentityTransformation())
    inv_transform = inv(transform)
    pmin, pmax = inv_transform.(extent(bounds))

    ax = range.(min.(pmin, pmax), max.(pmin, pmax), step=Δ)

    DoseGrid(ax...)
end

Base.size(pos::DoseGrid) = tuple(length.(pos.axes)...)
Base.length(pos::DoseGrid) = prod(length.(pos.axes))

Base.eachindex(pos::DoseGrid) = Base.OneTo(length(pos))
Base.CartesianIndices(pos::DoseGrid) = CartesianIndices(size(pos))

Base.getindex(pos::DoseGrid, i::Int) = pos[CartesianIndices(pos)[i]]

for op in (:+, :-, :*, :/)
    eval(quote
        function Base.$op(pos::DoseGrid, v::AbstractVector{<:Real})
            x = ($op).(getaxes(pos, 1), v[1])
            y = ($op).(getaxes(pos, 2), v[2])
            z = ($op).(getaxes(pos, 3), v[3])
            DoseGrid(x, y, z)
        end
    end)
end

#-- IO 

"""
    write_vtk(filename::String, pos::DoseGrid, data::Union{Vararg, Dict})

Save `DoseGrid` to the VTK Image data (vti) format.
"""
function write_vtk(filename::String, pos::DoseGrid, data::Vararg)
    vtk_grid(filename, pos.axes...) do vtkfile
        for (key, value) in data
            vtkfile[key, VTKPointData()] = value
        end
    end
end
write_vtk(filename::String, pos::DoseGrid, data::Dict) =  write_vtk(filename, pos, data...)

function get_nrrd_props(pos::AbstractDoseGrid)
    origin = tuple(pos[1, 1, 1]...)
    x, y, z = getaxes(pos)
    directions = (step(x), 0, 0), (0, step(y), 0), (0, 0, step(z))
    
    Dict("space origin"=>origin,
         "space directions"=>directions,
         "space"=>"left-posterior-superior",
         "kinds"=>("domain", "domain", "domain"))
end

"""
    write_nrrd(filename::String, pos::DoseGrid, data::AbstractArray)

Write `DoseGrid` to the NRRD data format.
"""
function write_nrrd(filename::String, pos::DoseGrid, data::AbstractArray)
    props = get_nrrd_props(pos::AbstractDoseGrid)
    FileIO.save(filename, reshape(data, size(pos)), props=props)
end

"""
    save(file::HDF5.H5DataStore, pos::DoseGrid)

Store `DoseGrid` data to an HDF5 file/group
"""
function save(file::HDF5.H5DataStore, pos::DoseGrid)
    for (i, ax) in enumerate("xyz")
        file["pos/$ax"] = collect(getaxes(pos, i))
    end
    nothing
end

"""
    load(::Type{DoseGrid}, file::HDF5.H5DataStore)

Load `DoseGrid` data from an HDF5 file/group
"""
function load(::Type{DoseGrid}, file::HDF5.H5DataStore)
    
    x = read(file["pos/x"])
    y = read(file["pos/y"])
    z = read(file["pos/z"])

    DoseGrid(x, y, z)
end

#--- DoseGridMasked ----------------------------------------------------------------------------------------------------
"""
    DoseGridMasked

Cartesian Dose Grid with a mask to reduce the number of dose points used.
"""
struct DoseGridMasked{TVec<:AbstractVector} <: AbstractDoseGrid
    axes::SVector{3, TVec}
    indices::Vector{CartesianIndex{3}}
    cells::Vector{Vector{Int}}
    function DoseGridMasked(x, y, z, indices, cells)
        ax = x, y, z
        if(!all(@. typeof(ax) <: StepRangeLen))
            ax = collect.(ax)
        end
        new{typeof(ax[1])}(SVector(ax...), indices, cells)
    end
end

"""
    DoseGridMasked(Δ, bounds::AbstractBounds, transform=I)

Construct a `DoseGridMasked` with spacing `Δ` within `bounds`.

The mask applies to all points within `bounds`. Can specify a coordinate system
transformation in `transform` to generate dose points in a different coordinate
system to the bounds.
"""
function DoseGridMasked(Δ, bounds::AbstractBounds, transform=IdentityTransformation())

    inv_transform = inv(transform)
    pmin, pmax = inv_transform.(extent(bounds))

    x, y, z = range.(min.(pmin, pmax), max.(pmin, pmax), step=Δ)

    grid_index = CartesianIndices( (length(x), length(y), length(z)) )
    indices = CartesianIndex{3}[]
    linear_index = zeros(Int, length(x), length(y), length(z))

    # Points
    for i in grid_index
        p = transform(SVector(x[i[1]], y[i[2]], z[i[3]]))
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

    DoseGridMasked(x, y, z, indices, cells)
end

for op in (:+, :-, :*, :/)
    eval(quote
        function Base.$op(pos::DoseGridMasked, v::AbstractVector{<:Real})
            x = ($op).(getaxes(pos, 1), v[1])
            y = ($op).(getaxes(pos, 2), v[2])
            z = ($op).(getaxes(pos, 3), v[3])
            DoseGridMasked(x, y, z, pos.indices, pos.cells)
        end
    end)
end

Base.length(pos::DoseGridMasked) = length(pos.indices)
Base.size(pos::DoseGridMasked) = (length(pos),)

gridsize(pos::DoseGridMasked) = tuple(length.(pos.axes)...)

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

# VTK

"""
    vtk_create_cell(cell)

Return the VTK cell type for a list of cell indices.
"""
function vtk_create_cell(cell)
    n = length(cell)
    n==5 && return MeshCell(VTKCellTypes.VTK_PYRAMID, cell)
    n==6 && return MeshCell(VTKCellTypes.VTK_WEDGE, cell)
    n==8 && return MeshCell(VTKCellTypes.VTK_VOXEL, cell)
end

"""
    write_vtk(filename::String, pos::DoseGridMasked, data::Union{Vararg, Dict})

Save `DoseGridMasked` to the VTK Unstructured Grid (vtu) format.
"""
function write_vtk(filename::String, pos::DoseGridMasked, data::Vararg)
    points = [pos[i] for i in eachindex(pos)]
    cells = vtk_create_cell.(pos.cells)

    vtk_grid(filename, points, cells) do vtkfile
        for (key, value) in data
            vtkfile[key, VTKPointData()] = value
        end
    end
end
write_vtk(filename::String, pos::DoseGridMasked, data::Dict) =  write_vtk(filename, pos, data...)

"""
    write_nrrd(filename::String, pos::DoseGridMasked, data::AbstractArray; fillvalue=NaN)

Write `DoseGridMasked` to the NRRD data format.

Masked points in `DoseGridMasked` are filled with `fillvalue`, defaults to `NaN`.
"""
function write_nrrd(filename::String, pos::DoseGridMasked, data::AbstractArray; fillvalue=NaN)

    props = get_nrrd_props(pos)

    datagrid = fill(fillvalue, gridsize(pos)...)
    for i in eachindex(data)
        I = pos.indices[i]
        datagrid[I] = data[i]
    end
    
    FileIO.save(filename, datagrid, props=props)
end

# HDF5

"""
    save(file::HDF5.H5DataStore, pos::DoseGridMasked)

Store `DoseGridMasked` data to an HDF5 file/group
"""
function save(file::HDF5.H5DataStore, pos::DoseGridMasked)
    
    for (i, ax) in enumerate("xyz")
        file["pos/$ax"] = collect(getaxes(pos, i))
    end

    file["pos/indices"] = hcat(collect.(Tuple.(pos.indices))...)

    file["cells/index"] = cumsum(vcat(1, length.(pos.cells)))
    file["cells/cells"] = vcat(pos.cells...)

    nothing
end

"""
    load(::Type{DoseGridMasked}, file::HDF5.H5DataStore)

Load `DoseGridMasked` data from an HDF5 file/group
"""
function load(::Type{DoseGridMasked}, file::HDF5.H5DataStore)
    
    x = read(file["pos/x"])
    y = read(file["pos/y"])
    z = read(file["pos/z"])

    indices = read(file["pos/indices"])
    indices = [CartesianIndex(col...) for col in eachcol(indices)]

    index = read(file["cells/index"])
    cells = read(file["cells/cells"])
    cells = [cells[index[i]:index[i+1]-1] for i in 1:length(index)-1]

    DoseGridMasked(x, y, z, indices, cells)
end
