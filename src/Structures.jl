export load_structure_from_ply, transform, transform!

"""
    load_structure_from_ply(filename)

Load a .ply file into a mesh

Uses MeshIO to load the mesh from the file. Unfortunately, MeshIO loads into
GeometryBasics.Mesh, not Meshes.SimpleMesh. The rest of the code converts the
mesh into a Meshes.SimpleMesh.
From https://github.com/JuliaIO/MeshIO.jl/issues/67#issuecomment-913353708
"""
function load_structure_from_ply(filename)
	mesh = FileIO.load(filename)
	points = [Float64.(Tuple(p)) for p in Set(mesh.position)]
	indices = Dict(p => i for (i, p) in enumerate(points))
	connectivities = map(mesh) do el
		Meshes.connect(Tuple(indices[Tuple(p)] for p in el))
	end
	Meshes.SimpleMesh(points, connectivities)
end

#--- Mesh Transformations -----------------------------------------------------
"""
    transform(mesh::SimpleMesh, trans)

Apply general transformation `trans` to `mesh``.
"""
function transform(mesh::SimpleMesh, trans)
    points = Point.(trans.(coordinates.(vertices(mesh))))
    SimpleMesh(points, topology(mesh))
end

"""
    transform!(mesh::SimpleMesh, trans)

Apply general transformation `trans` to `mesh`, modifying the original mesh.
"""
function transform!(mesh, trans)
    mesh.points .= Point.(trans.(coordinates.(vertices(mesh)))) 
end

#--- Mesh Intersection --------------------------------------------------------

"""
    intersect_mesh(s, mesh)

Find the intersection points of the `line` and the `mesh`.

Returns a list of intersection points, and an empty list if none present.
"""
function intersect_mesh(line::Geometry, mesh::Domain{Dim, T}) where {Dim, T}
    if(Threads.nthreads()==1)
        return intersect_mesh_single_threaded(line, mesh)
    else
        return intersect_mesh_multi_threaded(line, mesh)
    end
end

"""
    intersect_mesh_single_threaded(s, mesh)

Single-threaded version of `intersect_mesh`.
"""
function intersect_mesh_single_threaded(line::Geometry, mesh::Domain{Dim, T}) where {Dim, T}
    intersection_points = Point{Dim, T}[]
    intersection_cells = Geometry{Dim, T}[]
    for cell in mesh
        pt = intersect(line, cell)
        if(pt !== nothing)
            push!(intersection_points, pt)
            push!(intersection_cells, cell)
        end
    end
    intersection_points, intersection_cells
end

"""
    intersect_mesh_multi_threaded(s, mesh)

Multi-threaded version of `intersect_mesh`.
"""
function intersect_mesh_multi_threaded(line::Geometry, mesh::Domain{Dim, T}) where {Dim, T}
    intersection_points = Vector{Vector{Point{Dim, T}}}(undef, Threads.nthreads())
    intersection_cells = Vector{Vector{Geometry{Dim, T}}}(undef, Threads.nthreads())
    for i=1:Threads.nthreads()
        intersection_points[i] = Point{Dim, T}[]
        intersection_cells[i] = Geometry{Dim, T}[]
    end
    Threads.@threads for cell in mesh
        pt = intersect(line, cell)
        if(pt !== nothing)
            push!(intersection_points[Threads.threadid()], pt)
            push!(intersection_cells[Threads.threadid()], cell)
        end
    end
    vcat(intersection_points...), vcat(intersection_cells...)
end

#--- File IO ------------------------------------------------------------------

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


function isinside(testpoint::Point{3,T}, mesh::Mesh{3,T}) where T
    if !(eltype(mesh) <: Triangle)
      error("This function only works for surface meshes with triangles as elements.")
    end
    ex = testpoint .- extrema(mesh)
    direction = ex[argmax(norm.(ex))]
    r = Ray(testpoint, direction*2)
    
    intersects = false
    for elem in mesh
      if intersection(x -> x.type == NoIntersection ? false : true, r, elem)
        intersects = !intersects
      end
    end
    intersects
end
