#=
    Mesh Intersections

These tests load a sample mesh from an stl file, then compares the output of
`intersect_mesh`` to predefined intersection points.
=#

@testset "Meshes" begin
    mesh = load_structure_from_ply("test-data/test-mesh.stl")
    tol = 0.1

    test_intersection_point(p, p_truth, ε) = norm(p - p_truth) < ε

    function test_intersections(s, mesh, pts_truth, ε)
        pts = DoseCalculations.intersect_mesh(s, mesh)
        @test length(pts) == length(pts_truth)

        s = sortperm(getindex.(coordinates.(pts), 1))
        pts = pts[s]
        sT = sortperm(getindex.(coordinates.(pts_truth), 1))
        pts_truth = pts_truth[sT]

        @testset "$pI, $pIT" for (pI, pIT) in zip(pts, pts_truth)
            @test test_intersection_point(pI, pIT, ε)
        end
    end
    

    @testset "Mesh Intersections" begin

        @testset "Single Intersection" begin
            p1 = (-271.5, 30.2, -28.5)
            p2 = (152., 102., -52.)
            s = Segment(p1, p2)

            pts_truth = [Point(-166.7, 47.9, -34.3)]
            test_intersections(s, mesh, pts_truth, tol)
        end

        @testset "Double Intersection" begin
            p1 = (-299.9, -50.6, 139.3)
            p2 = (273.7, 113., -151.6)
            s = Segment(p1, p2)

            pts_truth = [Point(-134.77, -3.50, 55.56)
                         Point(159.98, 80.57, -93.93)]
                        
            test_intersections(s, mesh, pts_truth, tol)
        end

        @testset "No Intersection" begin
            p1 = (-299.9, -50.6, 139.3)
            p2 = (206., 100., 164.)
            s = Segment(p1, p2)

            pts = DoseCalculations.intersect_mesh(s, mesh)
            @test length(pts) == 0
        end

    end

    @testset "Partitioned Mesh Intersections" begin

        widths = DoseCalculations._boxwidth(mesh)./2

        parts = partition(mesh, BlockPartition(widths...))

        p1 = (-299.9, -50.6, 139.3)
        p2 = (273.7, 113., -151.6)
        s = Segment(p1, p2)
        pts_truth = DoseCalculations.intersect_mesh(s, mesh)

        test_intersections(s, parts, pts_truth, atol(Float64))
    end

    @testset "Closest Intersection" begin

        p1 = SVector(-299.9, -50.6, 139.3)
        p2 = SVector(273.7, 113., -151.6)

        pt_truth = SVector(-134.77, -3.50, 55.56)

        pt = closest_intersection(p1, p2, mesh)

        @test test_intersection_point(pt, pt_truth, tol)
        
    end
end
