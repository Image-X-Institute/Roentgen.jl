@testset "GantryPosition" begin
    ϕg, θb, SAD = rand(3)
    gantry = GantryPosition(ϕg, θb, SAD)

    a = [sin(ϕg), 0., cos(ϕg)]

    @testset "Constructor" begin
        @test gantry.gantry_angle ≈ ϕg
        @test gantry.collimator_angle ≈ θb
        @test gantry.source_axis_distance ≈ SAD
        @test gantry.central_beam_axis ≈ a
    end

    @testset "Methods" begin
        @test getϕg(gantry) ≈ ϕg
        @test getθb(gantry) ≈ θb
        @test getSAD(gantry) ≈ SAD
        @test DoseCalculations.beamaxis(gantry) ≈ a

        @test DoseCalculations.getposition(gantry) ≈ a*SAD

        fixed_to_bld(gantry) ≈ fixed_to_bld(ϕg, θb, SAD)
        bld_to_fixed(gantry) ≈ bld_to_fixed(ϕg, θb, SAD)
    end

end
