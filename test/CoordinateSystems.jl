@testset "CoordinateSystems" begin

    compare(v1, v2) = all(v1 .≈ v2)

    @testset "Fixed -> Gantry" begin
        pf = Vec(rand(3)...)
        @testset "ϕg = 0" begin
            ϕg = 0.
            pg = fixed_to_gantry(ϕg)(pf)
            @test compare(pf, pg)
            @test compare(pf, gantry_to_fixed(ϕg)(pg))
        end
        @testset "ϕg = 90" begin
            ϕg = 0.5*π
            pg = fixed_to_gantry(ϕg)(pf)
            @test compare(pf, Vec(pg[3], pg[2], -pg[1]))
            @test compare(pf, gantry_to_fixed(ϕg)(pg))
        end
    end

    @testset "Gantry -> BLD" begin
        pg = Vec(rand(3)...)
        SAD = rand()
        @testset "θb = 0" begin
            θb = 0.
            pb = gantry_to_bld(θb, SAD)(pg)
            @test compare(pg, pb + Vec(0, 0, SAD))
            @test compare(pg, bld_to_gantry(θb, SAD)(pb))
        end
        @testset "θb = 90" begin
            θb = 0.5*π
            pb = gantry_to_bld(θb, SAD)(pg)
            @test compare(pg, Vec(-pb[2], pb[1], pb[3] + SAD))
            @test compare(pg, bld_to_gantry(θb, SAD)(pb))
        end
    end

    @testset "Fixed -> BLD" begin
        pf = Vec(rand(3)...)
        SAD = rand()
        @testset "ϕg = 0, θb = 0" begin
            ϕg, θb = 0., 0.
            pb = fixed_to_bld(ϕg, θb, SAD)(pf)
            @test compare(pf, pb + Vec(0., 0., SAD))
            @test compare(pf, bld_to_fixed(ϕg, θb, SAD)(pb))
        end
        @testset "ϕg = 90, θb = 90" begin
            ϕg, θb = 0.5*π, 0.5*π
            pb = fixed_to_bld(ϕg, θb, SAD)(pf)
            @test compare(pf, Vec(pb[3] + SAD, pb[1], pb[2]))
            @test compare(pf, bld_to_fixed(ϕg, θb, SAD)(pb))
        end
        @testset "Inverse" begin
            ϕg, θb = 2*π*rand(), (2*rand()-1)*π
            pb = fixed_to_bld(ϕg, θb, SAD)(pf)
            @test all(pf .≈ bld_to_fixed(ϕg, θb, SAD)(pb))
        end
    end

end
|