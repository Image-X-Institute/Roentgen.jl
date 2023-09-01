#=
    Depth and SSD Calculation

These tests computes the Source-Surface Distance (SSD) and depth of two points
and compares their value to precomputed values. The two points are located
on-axis and off-axis. It also compares SSD "scaling": points that are along the
same ray line return the same SSD.

Implemented Surfaces:
    - ConstantSurface
    - PlaneSurface
    - MeshSurface (uses the same mesh and visual inspection as detailed in meshes.jl)
=#

@testset "External Surfaces" begin

    _test_mesh_path = "test-data/test-mesh.stl"

    function random_source(SAD)
        ϕ = 2π*rand()
        θ = π*rand()
        SAD*SVector(sin(ϕ)*cos(θ), cos(ϕ)*cos(θ), sin(θ))
    end

    random_position() = SVector((200*rand(3) .- 100)...)

    function test_surface(surf, pos, src, SSD_truth::T, depth_truth; atol=atol(T)) where T<:Real
        @test getSSD(surf, pos, src) ≈ SSD_truth atol=atol
        @test getdepth(surf, pos, src) ≈ depth_truth atol=atol
    
        λ = 1.05
        pos2 = src + λ*(pos - src)
        @test getSSD(surf, pos2, src) ≈ SSD_truth atol=atol
    end

    function test_surface(surf, pos, src, SSD_truth; args...)
        depth_truth = norm(pos-src)-SSD_truth
        test_surface(surf, pos, src, SSD_truth, depth_truth; args...)
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
        structure = load_structure_from_ply(_test_mesh_path)
        surf = MeshSurface(structure)

        @testset "Visual Inspection 1" begin
            src = SVector(0., 0., 1000.)
            pos = SVector(0., 0., 0.)
            test_surface(surf, pos, src, 902.4844182019233)
        end

        @testset "Visual Inspection 2" begin
            src = SVector(998.88987496197, 0., 47.10645070964268)
            pos = SVector(152., 102., -52.)

            test_surface(surf, pos, src, 846.4940339402172)
        end
    end

    @testset "Cylindrical Surface" begin
        mesh = load_structure_from_ply(_test_mesh_path)
        meshsurf = MeshSurface(mesh)
        surf = CylindricalSurface(mesh, 10., 13)
        write_vtk("surf", surf)
        @testset "Axis-Aligned, Center" begin
            src = SVector(0., 0., 1000.)
            pos = SVector(0., 0., 0.)
    
            test_surface(surf, pos, src, getSSD(meshsurf, pos, src); atol=1.)
        end
        
        @testset "Random position and source" begin
            src = SVector(998.88987496197, 0., 47.10645070964268)
            pos = SVector(52., 102., -52.)
            test_surface(surf, pos, src, getSSD(meshsurf, pos, src); atol=2.)
        end
    
        @testset "isinside" begin
            ϕ = range(0, 2π, length=7)
            y = range(-1, 1, length=5)
            rho = 1 .+ rand(length(ϕ), length(y))
            @. rho[end, :] = rho[1, :]
            pc = SVector(rand(3)...)
        
            surf = CylindricalSurface(ϕ, y, rho, pc)
        
            ϕᵢ = 2π*rand()
            yᵢ = (y[end]-y[1])*rand()+y[1]
        
            ρᵢ = Roentgen._interp(surf, ϕᵢ, yᵢ)
        
            from_cyl_coords(rho, ϕ, y) = SVector(rho*cos(ϕ), y, rho*sin(ϕ)) + pc

            ρ_inside = rand()*ρᵢ
            ρ_outside = rand()+ρᵢ
            y_inside = yᵢ
            y_outside_m = y[1] - rand()
            y_outside_p = y[end] + rand()

            @testset "Test positions" begin
                @test ρ_inside < ρᵢ
                @test ρ_outside > ρᵢ
                @test y[1]<=y_inside<y[end]
                @test y_outside_m < y[1]
                @test y_outside_p > y[end]                
            end
        
            @testset "Inside" begin
                @test Roentgen.isinside(surf, pc)
                @test Roentgen.isinside(surf, from_cyl_coords(ρ_inside, ϕᵢ, y_inside))
            end 
        
            @testset "Outside Radially" begin
                @test !Roentgen.isinside(surf, from_cyl_coords(ρ_outside, ϕᵢ, y_inside))
            end
        
            @testset "Outside Axially" begin
                @test !Roentgen.isinside(surf, from_cyl_coords(ρ_inside, ϕᵢ, y_outside_m))
                @test !Roentgen.isinside(surf, from_cyl_coords(ρ_inside, ϕᵢ, y_outside_p))
            end
        
            @testset "Outside Both" begin
                @test !Roentgen.isinside(surf, from_cyl_coords(ρ_outside, ϕᵢ, y_outside_m))
                @test !Roentgen.isinside(surf, from_cyl_coords(ρ_outside, ϕᵢ, y_outside_p))
            end
        end
    end

    @testset "Plane" begin
        n = SVector(rand(3)...)
        p = SVector(rand(3)...)
        plane = Roentgen.Plane(p, n)
    
        l1 = SVector(rand(3)...)
        l2 = SVector(rand(3)...)
        v = l2-l1
    
        @testset "Intersection" begin
            pI = Roentgen.intersection_point(plane, l1, l2)
    
            # On Plane
            @test dot(n, pI) ≈ dot(n, p)
    
            # On Line
            λ = (pI - l1)./v
            @test all(λ .≈ λ[1])
        end
    
        @testset "No Intersection" begin
            n = cross(SVector(rand(3)...), v)
            plane = Roentgen.Plane(p, n)
            @test Roentgen.intersection_point(plane, l1, l2) === nothing
        end
    end

    @testset "LinearSurface" begin

        function test_angle(surf, ϕg, SAD, SSDc)
        
            gantry = GantryPosition(ϕg, 0., SAD)
            src = Roentgen.getposition(gantry)
        
            T = RotY(ϕg)
        
            @testset "Along Central Axis" begin
                z = 20*rand().-10
                pos = T*SVector(0., 0., z)
                @test getSSD(surf, pos, src) ≈ SSDc
                @test getdepth(surf, pos, src) ≈ norm(pos-src)-SSDc
            end
        
            @testset "Off-Axis" begin
        
                R = 750.
                SSD = SSDc*√(1+(R/SAD)^2)
                
                θ = 2π*rand()
                z = 20*rand().-10
                α = (SAD-z)/SAD
        
                pos = T*SVector(α*R*sin(θ), α*R*cos(θ), z)
                @test getSSD(surf, pos, src) ≈ SSD
                @test getdepth(surf, pos, src) ≈ norm(pos-src)-SSD
            end
        end

        SAD = 1000.
        
        ϕg = [0., π/3., π, 4*π/3, 2π]
        r = @. 50*cos(ϕg) + 550.
        x = @. r*sin(ϕg)
        z = @. r*cos(ϕg)
        
        SSDc = SAD.-r
        
        p = SVector.(x, 0., z)
        n = normalize.(p)
        
        surf = LinearSurface(ϕg, p, n)

        @testset "$(rad2deg(ϕg[i])), $(SSDc[i])" for i in eachindex(ϕg, SSDc)
            test_angle(surf, ϕg[i], SAD, SSDc[i])
        end
        
    end
end
