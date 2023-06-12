@testset "Bixels" begin

    @testset "Constructors" begin
        x, y, wx, wy = rand(4)
        bixel = Bixel(x, y, wx, wy)
        @test bixel.position == [x, y]
        @test bixel.width == [wx, wy]

        x, y, w = rand(3)
        bixel = Bixel(x, y, w, w)
        @test bixel.position == [x, y]
        @test bixel.width == [w, w]

        x, w = rand(2), rand(2)
        bixel = Bixel(x, w)
        @test bixel.position == x
        @test bixel.width == w
        
    end

    @testset "Methods" begin
        pos = SVector(1., 2.)
        width = SVector(3., 4.)
        bixel = Bixel(pos, width)

        @testset "Center" begin
            @test getcenter(bixel) == pos
            @test getcenter(bixel, 1) == pos[1]
            @test bixel[1] == pos[1]
        end

        @testset "Width" begin
            @test getwidth(bixel) == width
            @test getwidth(bixel, 1) == width[1]
        end

        @testset "Edge" begin
            @test getedge(bixel) == pos - 0.5*width
            @test getedge(bixel, 1) == pos[1] - 0.5*width[1]
        end

        @test getarea(bixel) == width[1]*width[2]

        @testset "Subdivision" begin
            function test_subbixels(bixel, subbixels, n, subwidth)
                @test size(subbixels) == n

                @test all(getwidth.(subbixels) .≈ Ref(subwidth))
                @test sum(getarea.(subbixels)) ≈ getarea(bixel)

                # Test left-most edge of subbixels in both dimensions
                @test all(getedge.(subbixels[1, :], 1) .≈ getedge(bixel, 1))
                @test all(getedge.(subbixels[:, 1], 2) .≈ getedge(bixel, 2))

                # Test spacing between subbixels
                @test all(diff(getcenter.(subbixels, 1), dims=1) .≈ subwidth[1])
                @test all(diff(getcenter.(subbixels, 2), dims=2) .≈ subwidth[2])
            end

            n = 4, 3
            subwidth = getwidth(bixel)./n

            subbixels = subdivide(bixel, n...)
            test_subbixels(bixel, subbixels, n, subwidth)

            @test all(subdivide(bixel, subwidth...) .== subbixels)

        end
    end
end

@testset "Bixel Grids" begin
    # Default Constructor
    x = -12.:2.:8.
    y = -3.:3.:9.
    
    function test_bixelgrid(bixels, x, y)
        @test size(bixels) == (length(x)-1, length(y)-1)
        @test all(getwidth.(bixels) .== Ref([step(x), step(y)]))
        @test all(getedge.(bixels, 1) .== x[1:end-1])
        @test all(getedge.(bixels, 2) .== y[1:end-1]')
    end
    @testset "Default Constr." begin
        bixels = BixelGrid(x, y)
        test_bixelgrid(bixels, x, y)
    end

    # Snapped Range
    @testset "Snapped Range" begin
        Δx, Δy = 5., 2.
        xsnap, ysnap = DoseCalculations.snapped_range.((x, y), (Δx, Δy))
    
        bixels = BixelGrid(x, y, Δx, Δy)

        test_bixelgrid(bixels, xsnap, ysnap)
        Δ = 10*rand()
        @test BixelGrid(x, y, Δ, Δ) == BixelGrid(x, y, Δ)
    end
    
    # From Jaws

    function test_bixelgrid_jaws(bixels, jaws)
        left_edge = getedge.(bixels[1, :], 1)
        @test all(left_edge .== jaws.x[1])
        right_edge = @. getedge(bixels[end, :], 1)+getwidth(bixels[end, :], 1)
        @test all(right_edge .== jaws.x[2])

        bottom_edge = getedge.(bixels[:, 1], 2)
        @test all(bottom_edge .== jaws.y[1])
        top_edge = @. getedge(bixels[:, end], 2)+getwidth(bixels[:, end], 2)
        @test all(top_edge .== jaws.y[2])

        @test sum(getarea.(bixels)) == getarea(jaws)
    end

    @testset "From Jaws" begin
        Δx, Δy = 5., 2.
        jaws = Jaws(-13., 9., -1., 9.)
        bixels = BixelGrid(jaws, 5., 2.)

        test_bixelgrid_jaws(bixels, jaws)

        Δ = 10*rand()
        @test BixelGrid(x, y, Δ, Δ) == BixelGrid(x, y, Δ)
    end

    @testset "From MLC" begin
        function rand_mlc(y)
            n = length(y)-1
            xB = rand(n)
            xA = xB .+ rand(n)
            x = vcat(xB', xA')
            @. x = 20*x-10.
            
            MultiLeafCollimator(x, y)
        end

        y = -60:5.:60.
        mlc = rand_mlc(y)
        
        jaws = Jaws(-13., 18., -28., 16., )
        bixels = BixelGrid(mlc, jaws, 5.)
        
        test_bixelgrid_jaws(bixels, jaws)

        @test all(@. getwidth(bixels, 1)[2:end-1, :] == Δ)

        j1, j2 = locate.(Ref(mlc), gety(jaws))
        @test all(@. getedge(bixels, 2)[:, 2:end] == y[j1+1:j2]')
    end

end
