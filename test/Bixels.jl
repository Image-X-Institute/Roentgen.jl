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
        @time bixel = Bixel(x, w)
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
    end
end

@testset "Bixel Grid" begin

    Δ = SVector(1., 2.)

    x = -10.:Δ[1]:10.
    y = -12.:Δ[2]:10.
    grid = DoseCalculations.BixelGrid(x, y)

    nx = length(x)-1
    ny = length(y)-1

    @testset "Grid Size" begin
        @test size(grid) == (nx, ny)
        @test length(grid) == nx*ny
    end

    p = SVector(grid.xc[1], grid.yc[2])
    w = SVector(grid.Δx[1], grid.Δy[2])

    @testset "Position" begin
        @test getposition(grid, 1, 2) == p
        @test getwidth(grid, 1, 2) == w
    end

    b = Bixel(p, w)
    @testset "Indexing" begin
        grid[1, 2] == b # 2D Indexing
        grid[CartesianIndex(1, 2)] == b # Cartesian Indexing
        grid[LinearIndices(grid)[1, 2]] == b    # Linear Indexing
    end
end
