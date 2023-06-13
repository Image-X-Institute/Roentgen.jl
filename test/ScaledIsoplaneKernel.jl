#=
    (DISABLED) Tests for the ScaledIsoplaneKernel Dose Calculation Algorithm

See src/DoseCalculationAlgorithms/ScaledIsoplaneKernel.jl
=#
@testset "ScaledIsoplaneKernel" begin

    @testset "Sub-division" begin

        function test_subdivision(x, Δx, δxmax)
            x0, δx, nx = DoseCalculations.subdivide(x, Δx, δxmax)
            xsub = @. x0 + δx*(0:nx-1)
            
            @test δx <= δxmax
            @test x - 0.5*Δx ≈ xsub[1] - 0.5*step(xsub)
            @test x + 0.5*Δx ≈ xsub[end] + 0.5*step(xsub)

            x0, δx, nx
        end

        @testset "Odd" begin test_subdivision(0., 5., 1.) end
        @testset "Even" begin test_subdivision(0., 10., 1.) end
        @testset "Random" begin test_subdivision(0.48, 0.86, 0.78) end

        function test_subdivision_too_big(x, Δx, δxmax)
            x0, δx, nx = test_subdivision(x, Δx, δxmax)
            @test x0 ≈ x
            @test δx ≈ Δx
            @test nx == 1
        end

        @testset "Sub-division too big" begin test_subdivision_too_big(0., 1., 10.) end
        @testset "Sub-division equal" begin test_subdivision_too_big(0., 2., 2.) end

    end

    @testset "Sub-cycling" begin

        jaws = Jaws(10.)
        
        pos = [SVector(0., 0., -1000.)]
        surf = PlaneSurface(800.)

        calc = ScaledIsoplaneKernel("src/DoseCalculationAlgorithms/sample-kernel-data/6x/", 280.)

        function get_dose(Δx, Δy)
            bixels = DoseCalculations.BixelGrid(jaws, Δx, Δy)
            sum(dose_fluence_matrix(pos, vec(bixels), surf, calc))
        end

        ref_dose = get_dose(1., 1.)

        gridsizes = [(5., 5.),
                     (5., 1.),
                     (1., 5.)]

        @testset for (Δx, Δy) in gridsizes
            @test ref_dose ≈ get_dose(Δx, Δy)
        end
    end

end