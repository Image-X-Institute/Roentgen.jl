@testset "Fluence" begin
    @testset "Overlap" begin
        x = 10*rand()
        w = 2*rand()

        @testset "Bixel within range" begin
            @test DoseCalculations.overlap(x, w, x-2*w, x+2*w) ≈ 1.
        end
        @testset "Bixel outside range" begin
            @test DoseCalculations.overlap(x, w, x+2*w, x+3*w) ≈ 0.
            @test DoseCalculations.overlap(x, w, x-3*w, x-2*w) ≈ 0.
        end
        @testset "Bixel partially within range" begin
            @test DoseCalculations.overlap(x, w, x-2*w, x+0.1*w) ≈ 0.6
            @test DoseCalculations.overlap(x, w, x+0.1*w, x+2*w) ≈ 0.4
        end
        @testset "Bixel fully within range" begin
            @test DoseCalculations.overlap(x, w, x-0.1*w, x+0.1*w) ≈ 0.2
        end
    end
end