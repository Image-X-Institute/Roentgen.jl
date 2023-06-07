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

    @testset "Fluence from Rectangle" begin
        pos, width = [0.5, 0.], [1., 2.]
        bixel = Bixel(pos, width)

        # Bixel fully inside
        @test DoseCalculations.fluence_from_rectangle(bixel, [-0.2, 1.2], [-1.2, 1.2]) ≈ 1.
        
        # x half inside
        x = rand()
        @test DoseCalculations.fluence_from_rectangle(bixel, [-0.2, x], [-1.2, 1.2]) ≈ x/getwidth(bixel, 1)
        x = rand()
        @test DoseCalculations.fluence_from_rectangle(bixel, [x, 1.2], [-1.2, 1.2]) ≈ (1-x)/getwidth(bixel, 1)
        
        # y half inside
        y = 2*rand()-1
        @test DoseCalculations.fluence_from_rectangle(bixel, [-0.2, 1.2], [-1.2, y]) ≈ (y+1)/getwidth(bixel, 2)
        y = 2*rand()-1
        @test DoseCalculations.fluence_from_rectangle(bixel, [-0.2, 1.2], [y, 1.2]) ≈ (1-y)/getwidth(bixel, 2)

        # x and y inside
        x = sort(rand(2))
        y = sort(2*rand(2).-1)
        @test DoseCalculations.fluence_from_rectangle(bixel, x, y) ≈ (x[2]-x[1])*(y[2]-y[1])/getarea(bixel)
    end
end
