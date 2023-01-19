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

    function test_surface(surf, pos, src, SSD_truth, depth_truth)
        @test getSSD(surf, pos, src) ≈ SSD_truth
        @test getdepth(surf, pos, src) ≈ depth_truth
    
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

        test_surface(surf, pos, src, SSD, d)
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

        @testset "3-4-5 Triangle" test_surface(surf, pos, src, SSD, d)

        # Rotationally Invariant
        @testset "Rotational Invariance" begin
            T = RotXYZ(RotXYZ(2π*rand(3)...))            
            test_surface(surf, T*pos, T*src, SSD, d)
        end
    end

    @testset "MeshSurface" begin
        structure = load_structure_from_ply("test/test_mesh.stl")
        surf = MeshSurface(structure)

        # Test 1 - Visually inspected for accuracy

        @testset "Visual Inspection 1" begin
            src = SVector(0., 0., 1000.)
            pos = SVector(0., 0., 0.)
            test_surface(surf, pos, src, 884.0906064830797, 115.90939351692032)
        end

        # Test 2 - Visually inspected for accuracy

        @testset "Visual Inspection 2" begin
            src = SVector(-335, 0., 942)
            pos = SVector(30., 20., 10.)

            test_surface(surf, pos, src, 875.0481662974585, 126.075702162384)
        end 
    end

    @testset "VariablePlaneSurface" begin
        SAD = 1000.

        @testset "Constant Variable Surface" begin
            SSD₀ = 800.

            ϕ = deg2rad.(-180:180)
            distance = fill(SSD₀, length(ϕ))

            surf = VariablePlaneSurface(ϕ, distance)
            planesurf = PlaneSurface(SSD₀)

            src = random_source(SAD)
            pos = random_position()

            SSD = getSSD(planesurf, pos, src)
            d = getdepth(planesurf, pos, src)

            test_surface(surf, pos, src, SSD, d)
        end

        @testset "Variable Surface" begin

            ϕ = deg2rad.(-180:180)
            distance = 750. .+ 200*rand(length(ϕ))

            surf = VariablePlaneSurface(ϕ, distance)

            src = random_source(SAD)
            pos = random_position()

            SSD₀ = DoseCalculations.interpolate(surf, src)
            planesurf = PlaneSurface(SSD₀)

            SSD = getSSD(planesurf, pos, src)
            d = getdepth(planesurf, pos, src)

            test_surface(surf, pos, src, SSD, d)
        end
    end

end

