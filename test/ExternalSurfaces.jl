using Meshes, LinearAlgebra, StaticArrays

#=
    Depth and SSD Calculation

These tests computes the Source-Surface Distance (SSD) and depth of two points
and compares their value to precomputed values. The two points are located
on-axis and off-axis. It also compares SSD "scaling": points that are along the
same ray line return the same SSD.

Implemented Surfaces:
    - ConstantSurface
    - PlaneSurface
    - MeshSurface (uses the same mesh and visual inspection as detailed in
      test/meshes.jl)
=#

function test_external_surfaces(surf, test_points, ε=1e-3)
    @testset "On Axis" begin
        @test ≈(getSSD(surf, test_points[1]["pos"]), test_points[1]["SSD"], atol=ε)
        @test ≈(getdepth(surf, test_points[1]["pos"]), test_points[1]["depth"], atol=ε)
    end
    @testset "Off Axis" begin
        @test ≈(getSSD(surf, test_points[2]["pos"]), test_points[2]["SSD"], atol=ε)
        @test ≈(getdepth(surf, test_points[2]["pos"]), test_points[2]["depth"], atol=ε)
    end

    @testset "SSD Scaling" begin
        x, y, z = 19., -17., 1000. # Positions and scaling chosen at random
        α = 2.31
        p₁ = SVector(x, y, z)
        p₂ = SVector(α*x, α*y, α*z)

        @test getSSD(surf, p₁) ≈ getSSD(surf, p₂)
        @test getdepth(surf, p₂) - getdepth(surf, p₁) ≈ norm(p₁ - p₂)
    end
end

@testset "External Surfaces" begin

    @testset "ConstantSurface" begin
        surf = ConstantSurface(800.)

        test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>800., "depth"=>200.),
                       Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>800., "depth"=>558.9930095478785)]
        # Off axis x and y positions chosen at random

        test_external_surfaces(surf, test_points)
    end

    @testset "PlaneSurface" begin
        surf = PlaneSurface(800.)

        test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>800., "depth"=>200.),
                       Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>802.3574964120315, "depth"=>556.635513135847)]

        test_external_surfaces(surf, test_points)
    end

    @testset "MeshSurface" begin
        structure = load_structure_from_ply("test_mesh.stl")
        for (i, point) in enumerate(vertices(structure))
            structure.points[i] =  convert(Point{3, Float64}, point + Vec(0., 0., 1000.))
        end
        surf = MeshSurface(structure)

        test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>884.0906064830797, "depth"=>115.90939351692032),
                       Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>929.7027754972272, "depth"=>429.2902340506513)]

        test_external_surfaces(surf, test_points)
    end

    @testset "IsoplaneSurface" begin
        structure = load_structure_from_ply("test_mesh.stl")
        SAD = 1000.
        mesh_surf = MeshSurface(transform(structure, fixed_to_bld(0., 0., SAD)))

        x = -40:1.:40
        y = -70:1.:70
        surf = IsoplaneSurface(x, y, 1000.)

        compute_SSD!(surf, mesh_surf)

        test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>884.0906064830797, "depth"=>115.90939351692032),
                       Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>929.7027754972272, "depth"=>429.2902340506513)]

        test_external_surfaces(surf, test_points)
    end
end

