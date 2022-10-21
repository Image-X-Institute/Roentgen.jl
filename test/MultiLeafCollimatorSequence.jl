@testset "MultiLeafCollimatorSequence" begin 

    function test_constructor(mlc, x, y)
        @test mlc.positions == x
        @test mlc.edges == y
        @test mlc.nleaves == length(y)-1
        @test 2 == size(x, 1)
        @test mlc.nleaves == size(x, 2)
        @test mlc.napertures == size(x, 3)
    end

    function test_methods(mlc)
        @test length(mlc) == mlc.napertures
        @test size(mlc) == (mlc.napertures,)
        @test firstindex(mlc) == 1
        @test lastindex(mlc) == mlc.napertures
        @test eachindex(mlc) == Base.OneTo(mlc.napertures)
    end

    function test_indexing(mlc)
        i = 3
        @test mlc[i] == MultiLeafCollimator(mlc.positions[:, :, i], mlc.edges)

        i = 2:6
        x = mlc.positions[:, :, i]
        @test mlc[i] == MultiLeafCollimatorSequence(x, mlc.edges)
    end

    function test_operations(mlc, Δ)
        x = mlc.positions
        y = mlc.edges
        
        @test mlc + Δ == MultiLeafCollimatorSequence(x.+Δ[1], y.+Δ[2])
        @test mlc - Δ == MultiLeafCollimatorSequence(x.-Δ[1], y.-Δ[2])
        @test mlc / Δ == MultiLeafCollimatorSequence(x./Δ[1], y./Δ[2])
        @test mlc * Δ == MultiLeafCollimatorSequence(x.*Δ[1], y.*Δ[2])
    end

    function test_io_hdf5(mlc)

        filename = "tmp.hdf5"
        h5open(filename, "w") do file
            save(file, mlc)
        end

        file = h5open(filename, "r")
        close(file)

        @testset "Write" begin
            h5open(filename, "r") do file
                @test haskey(file, "mlc")
                @test haskey(file, "mlc/edges")
                @test haskey(file, "mlc/positions")

                @test read(file["mlc/edges"]) == getedges(mlc)
                @test read(file["mlc/positions"]) == getpositions(mlc)
            end
        end

        @testset "Read" begin
            mlc2 = h5open(filename, "r") do file
                load(MultiLeafCollimatorSequence, file)
            end
            @test mlc2 == mlc
        end

        rm(filename)

    end

    y = -30.:5:30.
    n = 10
    @testset "With $(typeof(yᵢ)) edges" for yᵢ in [y, collect(y)]
        x = rand(2, length(yᵢ)-1, n)
        mlc = MultiLeafCollimatorSequence(x, yᵢ)
        
        @testset "Constructor" test_constructor(mlc, x, yᵢ)

        @testset "Constructor Errors" test_constructor(mlc, x, yᵢ)

        @testset "Basic Methods" test_methods(mlc)
        @testset "Indexing" test_indexing(mlc)
        
        @testset "Operate with Vector" test_operations(mlc, rand(2))
        @testset "Operate with Tuple" test_operations(mlc, tuple(rand(2)...))

        @testset "HDF5 - IO" test_io_hdf5(mlc)
    end

    y = -30.:5:30.
    @testset "Assertion Errors" begin
        x = rand(2, length(y)-10, 3)
        @test_throws AssertionError("Length of positions and edges do not match") MultiLeafCollimatorSequence(x, y)

        x = rand(1, length(y)-1, 4)
        @test_throws AssertionError("Number of leaf positions per track != 2") MultiLeafCollimatorSequence(x, y)
    end

end
