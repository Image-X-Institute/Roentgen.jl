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
    pts = DoseCalculations.intersect_mesh(s, mesh)
    @test length(pts) == length(pts_truth)
    @test all(test_intersection_point.(pts, pts_truth, ε))
end

@testset "Mesh Intersections" begin
    mesh = load_structure_from_ply("test_mesh.stl")
    ε = 0.1

    @testset "Single Intersection" begin
        s = Segment((63., -15., 2.), (0., 0., 1000.))
        pts_truth = [Point([57.3, -13.6, 92.6])]
        test_intersections(s, mesh, pts_truth, ε)
    end

    @testset "Double Intersection" begin
        s = Segment((-73., 7., -191.), (0., 0., 1000.))

        pts_truth = [Point([-67.0, 6.4, -92.4]),
                     Point([-54.1, 5.2, 116.9])]
                    
        test_intersections(s, mesh, pts_truth, ε)
    end

    @testset "No Intersection" begin
        s = Segment((-189., 53.8, -136.3), (0., 0., 1000.))
        pts = DoseCalculations.intersect_mesh(s, mesh)
        @test length(pts) == 0
    end

end
