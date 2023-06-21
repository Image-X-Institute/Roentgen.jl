#=
    Mesh Intersections

These tests load a sample mesh from an stl file, then compares the output of
`intersect_mesh`` to predefined intersection points.

The mesh consists of an elliptical cylinder of radii 200mm and 300mm, which has
been arbitrarily rotated.

The starting positions of intersecting lines (segments) where chosen randomly
to be either within or outside the mesh, depending on the test. The line ends
at the source position (0, 0, 1000).

The true intersection points were obtained by visually inspecting intersecting
lines in a 3D visualisation program (there is a Paraview state file in `meshes/`
which displays the mesh, all lines and all intersecting points).

It tests three types of "intersections":
    - Only one intersection (point inside cylinder)
    - Double intersection (cylinder between source and point)
    - No intersection.
=#

test_intersection_point(p, p_truth, ε) = norm(p - p_truth) < ε

function test_intersections(s, mesh, pts_truth, ε)
    pts = DoseCalculations.intersect_mesh_single_threaded(s, mesh)
    @test length(pts) == length(pts_truth)
    @testset "$pI, $pIT" for (pI, pIT) in zip(pts, pts_truth)
        @test test_intersection_point(pI, pIT, ε)
    end
end

@testset "Mesh Intersections" begin
    mesh = load_structure_from_ply("test-data/test-mesh.stl")
    ε = 0.1

    @testset "Single Intersection" begin
        p1 = (-271.5, 30.2, -28.5)
        p2 = (152., 102., -52.)
        s = Segment(p1, p2)

        pts_truth = [Point(-166.7, 47.9, -34.3)]
        test_intersections(s, mesh, pts_truth, ε)
    end

    @testset "Double Intersection" begin
        p1 = (-299.9, -50.6, 139.3)
        p2 = (273.7, 113., -151.6)
        s = Segment(p1, p2)

        DoseCalculations.intersect_mesh_single_threaded(s, mesh)

        pts_truth = [Point(-134.77, -3.50, 55.56)
                     Point(159.98, 80.57, -93.93)]
                    
        test_intersections(s, mesh, pts_truth, ε)
    end

    @testset "No Intersection" begin
        p1 = (-299.9, -50.6, 139.3)
        p2 = (206., 100., 164.)
        s = Segment(p1, p2)

        pts = DoseCalculations.intersect_mesh(s, mesh)
        @test length(pts) == 0
    end

end
