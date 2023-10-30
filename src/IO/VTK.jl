
#--- DoseGrid ----------------------------------------------------------------------------------------------------------

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

#--- DoseGridMasked ---------------------------------------------------------------------------------------------------------

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
