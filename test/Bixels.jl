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

@testset "Bixel Grid" begin
    # Default Constructor
    
    function test_bixelgrid(bixels, x, y)

        @test size(bixels) == (length(x)-1, length(y)-1)

        # x
        @test all(getwidth.(bixels, 1) .== diff(x))
        @test all(getwidth.(bixels, 1) .> 0.)
        @test all(getedge.(bixels, 1) .== x[1:end-1])

        # y
        @test all(getwidth.(bixels, 2) .== diff(y)')
        @test all(getwidth.(bixels, 2) .> 0.)
        @test all(getedge.(bixels, 2) .== y[1:end-1]')
    end

    test_bixelgrid_methods(bixels) = test_bixelgrid(bixels, getaxes(bixels)...)

    @testset "Default Constr." begin
        x = -12.:2.:8.
        y = -3.:3.:9.
        bixels = BixelGrid(x, y)
        test_bixelgrid(bixels, x, y)
    end

    # Snapped Range
    @testset "Snapped Range" begin
        x = -12.:2.:8.
        y = -3.:3.:9.
        Δx, Δy = 5., 2.
        xsnap, ysnap = Roentgen.snapped_range.((x, y), (Δx, Δy))
    
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

        @testset "Constructor" test_bixelgrid_jaws(bixels, jaws)
        @testset "Methods" test_bixelgrid_methods(bixels)

        Δ = 10*rand()
        @test BixelGrid(jaws, Δ, Δ) == BixelGrid(jaws, Δ)
    end

end

@testset "Bixels from BLD" begin
    x = [-3. -20. 3. -5. 11. -15.
          2.  20. 3.  7. 13. -8.]
    y = -15.:5.:15.
    mlc = MultiLeafCollimator(x, y)
    jaws = Jaws(-13., 18., -12., 14., )
    
    test_data_path = "test-data/bixels_from_aperture.jld2"

    bixels = bixels_from_bld(mlc, jaws; Δx=5., snap_to_aperture=true)
    @test bixels == JLD2.load(test_data_path, "snap=true")

    bixels = bixels_from_bld(mlc, jaws; Δx=5., snap_to_aperture=false)
    @test bixels == JLD2.load(test_data_path, "snap=false")

    @testset "Issue #87" begin
        mlc, jaws = load("test-data/issue87.jld2", "mlc", "jaws")
        bixels = bixels_from_bld(mlc, jaws; snap_to_aperture=true)

        @test all(@. getwidth(bixels, 1)>0)
        @test all(@. getwidth(bixels, 2)>0)
    end
end
