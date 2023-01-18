using Meshes, LinearAlgebra, StaticArrays, Rotations

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

@testset "External Surfaces" begin

    function random_source(SAD)
        ϕ = 2π*rand()
        θ = π*rand()
        SAD*SVector(sin(ϕ)*cos(θ), cos(ϕ)*cos(θ), sin(θ))
    end

    function random_position()
        x, y, z = 200*rand(3) .- 100
        SVector(x, y, z)
    end



    function SSD_scale_invariance(surf, src, pos, λ=2.31)
        pos2 = src + λ*(pos - src)
        getSSD(surf, src, pos2) ≈ getSSD(surf, src, pos)
    end

    SAD = 1000.
    @test norm(random_source(SAD)) ≈ SAD

    @testset "ConstantSurface" begin
        SSD = 800.
        surf = ConstantSurface(SSD)

        src = random_source(SAD)
        pos = random_position()

        d = norm(pos - src) - SSD

        @test getSSD(surf, src, pos) == SSD
        @test getdepth(surf, src, pos) == d
        @test SSD_scale_invariance(surf, src, pos)

        # Rotationally Invariant
        T = RotY(2π*rand())
        @test getSSD(surf, T*src, T*pos) == SSD
        @test getdepth(surf, T*src, T*pos) == d
        @test SSD_scale_invariance(surf, T*src, T*pos)
    end

    # @testset "PlaneSurface" begin
    #     surf = PlaneSurface(800.)

    #     test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>800., "depth"=>200.),
    #                    Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>802.3574964120315, "depth"=>556.635513135847)]

    #     test_external_surfaces(surf, test_points)
    # end

    # @testset "MeshSurface" begin
    #     structure = load_structure_from_ply("test_mesh.stl")
    #     for (i, point) in enumerate(vertices(structure))
    #         structure.points[i] =  convert(Point{3, Float64}, point + Vec(0., 0., 1000.))
    #     end
    #     surf = MeshSurface(structure)

    #     test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>884.0906064830797, "depth"=>115.90939351692032),
    #                    Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>929.7027754972272, "depth"=>429.2902340506513)]

    #     test_external_surfaces(surf, test_points)
    # end

    # @testset "IsoplaneSurface" begin
    #     structure = load_structure_from_ply("test_mesh.stl")
    #     SAD = 1000.
    #     mesh_surf = MeshSurface(transform(structure, fixed_to_bld(0., 0., SAD)))

    #     x = -40:1.:40
    #     y = -70:1.:70
    #     surf = IsoplaneSurface(x, y, 1000.)

    #     compute_SSD!(surf, mesh_surf)

    #     test_points = [Dict("pos"=>SVector(0., 0., 1000.), "SSD"=>884.0906064830797, "depth"=>115.90939351692032),
    #                    Dict("pos"=>SVector(-54., 89., 1355.), "SSD"=>929.7027754972272, "depth"=>429.2902340506513)]

    #     test_external_surfaces(surf, test_points)
    # end
end

