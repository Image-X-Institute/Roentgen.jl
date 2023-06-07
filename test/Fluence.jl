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

    @testset "Fluence from Beam Limiting Devices" begin

        @testset "MultiLeafCollimator" begin

            mlcx = [-10. -5. -8.
                      5. 12.  6.]
            mlcy = -10.:5.:5
            
            mlc = MultiLeafCollimator(mlcx, mlcy)
            
            # Inside aperture
            bixel = Bixel(1., 0.1, 2., 3.)
            @test fluence(bixel, mlc) ≈ 1.
            
            # Outside aperture
            bixel = Bixel(-10., -3., 5., 3.)
            @test fluence(bixel, mlc) ≈ 0.
            
            # Overlapping aperture
            bixel = Bixel(5., -3., 6., 8.)
            
            Δx = [mlcx[2, 3] - getedge(bixel, 1)
                  getwidth(bixel, 1)
                  mlcx[2, 1] - getedge(bixel, 1)]
            
            Δy = [getedge(bixel, 2)+getwidth(bixel, 2)-mlcy[3],
                  mlcy[3]-mlcy[2],
                  mlcy[2]-getedge(bixel, 2)]
            @test fluence(bixel, mlc) ≈ sum(Δx.*Δy)/getarea(bixel)
            
            # Suppying an index            
            bixel = Bixel(6., 2.5, 5., 5.)
            @test fluence(bixel, mlc) ≈ fluence(bixel, 3, getpositions(mlc))

        end

        @testset "Jaws" begin
            
            bixel = Bixel(0.5, 0., 1., 2.)
            x = sort(rand(2))
            y = sort(2*rand(2).-1)
            jaws = Jaws(x, y)

            @test fluence(bixel, jaws) ≈ (x[2]-x[1])*(y[2]-y[1])/getarea(bixel)
        end
    end
end

using DoseCalculations
using Plots
using BeamLimitingDevicePlots
using Test

begin
    p = plot()
    plot_bld!(p, bixel)
    plot_bld!(p, mlc; invert=true)
end

using BenchmarkTools
@btime fluence($bixel, $mlc);
