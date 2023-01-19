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

    random_position() = SVector((200*rand(3) .- 100)...)

    function test_SSD_calculation(surf, pos, src, SSD_truth)
        @test getSSD(surf, pos, src) ≈ SSD_truth
    
        λ = 2*rand()
        pos2 = src + λ*(pos - src)
        @test getSSD(surf, pos2, src) ≈ SSD_truth
    end
    SAD = 1000.
    @test norm(random_source(SAD)) ≈ SAD

    @testset "ConstantSurface" begin
        SAD = 1000.
        SSD = 800.
        surf = ConstantSurface(SSD)

        src = random_source(SAD)
        pos = random_position()

        d = norm(pos - src) - SSD

        test_SSD_calculation(surf, pos, src, SSD)
        @test getdepth(surf, pos, src) ≈ d
    end

    @testset "PlaneSurface" begin
        SAD = 1000.

        # 3-4-5 Triangle
        SSD₀ = 400. # Central axis source-surface distance
        ρ = 300.

        surf = PlaneSurface(SSD₀)

        SSD = √(SSD₀^2 + ρ^2)

        ϕ = rand()*2*π
        z = 20*rand()-10
        ρ′ = ρ*(SAD-z)/SSD₀

        src = SAD*SVector(0., 0., 1.)
        pos = SVector(ρ′*cos(ϕ), ρ′*sin(ϕ), z)
        
        d = norm(pos - src) - SSD

        test_SSD_calculation(surf, pos, src, SSD)
        @test getdepth(surf, pos, src) ≈ d

        # Rotationally Invariant
        T = RotXYZ(RotXYZ(2π*rand(3)...))
        test_SSD_calculation(surf, T*pos, T*src, SSD)
        @test getdepth(surf, T*pos, T*src) ≈ d
    end

    @testset "MeshSurface" begin
        structure = load_structure_from_ply("test/test_mesh.stl")
        surf = MeshSurface(structure)

        # Test 1 - Visually inspected for accuracy

        src = SVector(0., 0., 1000.)
        pos = SVector(0., 0., 0.)

        test_SSD_calculation(surf, pos, src, 884.0906064830797)
        @test getdepth(surf, pos, src) ≈ 115.90939351692032

        # Test 2 - Visually inspected for accuracy

        src = SVector(-335, 0., 942)
        pos = SVector(30., 20., 10.)
        test_SSD_calculation(surf, pos, src, 875.0481662974585)
        @test getdepth(surf, pos, src) ≈ 126.075702162384
    end

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
