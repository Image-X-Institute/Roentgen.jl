
#= Dose Grids =========================================================================================================#

#--- DoseGrid -----------------------------------------------------------------

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
write_vtk(filename::String, pos::DoseGrid, data::Dict) = write_vtk(filename, pos, data...)

#--- DoseGridMasked -----------------------------------------------------------

"""
    write_vtk(filename::String, pos::DoseGridMasked, data::Union{Vararg, Dict})

Save `DoseGridMasked` to the VTK Unstructured Grid (vtu) format.
"""
function write_vtk(filename::String, pos::DoseGridMasked, data::Vararg)
    points = collect(pos)
    
    cells = _vtk_create_cell.(_get_cells(pos))

    vtk_grid(filename, points, cells) do vtkfile
        for (key, value) in data
            vtkfile[key, VTKPointData()] = value
        end
    end
end
write_vtk(filename::String, pos::DoseGridMasked, data::Dict) =  write_vtk(filename, pos, data...)

"""
    _get_cells(pos)

Returns a list of cell connectivities
"""
function _get_cells(pos)
    neighbours = vec([CartesianIndex(i,j,k) for i=0:1, j=0:1, k=0:1])[2:end]

    indices = CartesianIndices(pos)
    linear_index = zeros(gridsize(pos))

    for i in eachindex(pos)
        linear_index[indices[i]] = i
    end

    cells = Vector{Int}[]
    for i in eachindex(indices)
        cell = [i]
    
        for n in neighbours
            index = indices[i] + n
            if(checkbounds(Bool, linear_index, index) && linear_index[index]!=0)
                push!(cell, linear_index[index])
            end
        end
        if(length(cell)==8)
            push!(cells, cell)
        end
    end
    cells
end

"""
    _vtk_create_cell(cell)

Return the VTK cell type for a list of cell indices.
"""
function _vtk_create_cell(cell)
    n = length(cell)
    n==5 && return MeshCell(VTKCellTypes.VTK_PYRAMID, cell)
    n==6 && return MeshCell(VTKCellTypes.VTK_WEDGE, cell)
    n==8 && return MeshCell(VTKCellTypes.VTK_VOXEL, cell)
end

#= External Surfaces ==================================================================================================#

#--- MeshSurface --------------------------------------------------------------

write_vtk(filename::String, surf::MeshSurface) = write_vtk(filename, surf.mesh)

#--- CylindricalSurface -------------------------------------------------------

function write_vtk(filename::String, surf::CylindricalSurface)
    ϕ = surf.ϕ
    y = surf.y
    rho = Interpolations.coefficients(surf.rho)
    @. rho[isinf(rho)] = NaN
    xc, yc, zc = surf.center

    x = @. xc + rho*sin(ϕ)
    y = @. yc + surf.y
    z = @. zc + rho*cos(ϕ)

    xg = reshape(x, size(x)..., 1)
    yg = ones(size(xg)).*y'
    zg = reshape(z, size(z)..., 1)

    vtk = vtk_grid(filename, xg, yg, zg)
    vtk_save(vtk)
end

#= Structures =========================================================================================================#

"""
    write_vtk(filename, mesh::SimpleMesh)

Save a `SimpleMesh` to a VTK (.vtu) file.
"""
function write_vtk(filename, mesh::SimpleMesh)
    id = [indices(e) for e in elements(topology(mesh))]
    vtkfile = vtk_grid(filename,
                       coordinates.(vertices(mesh)), 
                       MeshCell.(Ref(VTKCellTypes.VTK_TRIANGLE), id))
    vtk_save(vtkfile)
end

#= Beamlets ===========================================================================================================#

"""
    _vtk_beamlet_points(beamlet, L)

Computes the vertices of the pyramid covered by a beamlet.

See `write_vtk(filename, beamlet::Beamlet, L=2)` for further details
"""
function _vtk_beamlet_points(beamlet, L)
    s = Roentgen.source_position(beamlet)
    D = L*Roentgen.source_axis_distance(beamlet)
    ax, ay, az = Roentgen.beamlet_axes(beamlet)
    wx, wy = L.*Roentgen.halfwidth(beamlet)

    (s,) .+ [D*az+wx*ax+wy*ay, D*az-wx*ax+wy*ay, D*az-wx*ax-wy*ay, D*az+wx*ax-wy*ay, zeros(3)]
end

"""
    write_vtk(filename, mesh::Beamlet, L=2)

Save a `Beamlet` to a VTK (.vtu) file.

The beamlet is drawn from the source position to a length determined by parameter
`L`. This is the length of the beamlet in units of source-axis distance.
"""
function write_vtk(filename, beamlet::Beamlet, L=2)
    points = _vtk_beamlet_points(beamlet, L)
    cells = [MeshCell(VTKCellTypes.VTK_PYRAMID, collect(1:5))]
    vtk = vtk_grid(filename, points, cells)
    vtk_save(vtk)
end

"""
    write_vtk(filename, mesh::Vector{<:AbstractBeamlet}, L=2)

Save a vector of `Beamlet` to a VTK (.vtu) file.

See `write_vtk(filename, beamlet::Beamlet, L=2)` for further details
"""
function write_vtk(filename, beamlets::Vector{<:AbstractBeamlet}, L=2)
    n = length(beamlets)

    points = SVector{3, Float64}[]
    for i in eachindex(beamlets)
        append!(points, _vtk_beamlet_points(beamlets[i], L))
    end

    indices =  reshape(1:n*5, 5, n)
    cells = MeshCell.((VTKCellTypes.VTK_PYRAMID,), eachcol(indices))

    vtk = vtk_grid(filename, points, cells)
    vtk_save(vtk)
end
