@testset "MultiLeafCollimator" begin 

    function test_operations(mlc, Δ)
        x = mlc.positions
        y = mlc.edges
        
        @test mlc + Δ == MultiLeafCollimator(x.+Δ[1], y.+Δ[2])
        @test mlc - Δ == MultiLeafCollimator(x.-Δ[1], y.-Δ[2])
        @test mlc / Δ == MultiLeafCollimator(x./Δ[1], y./Δ[2])
        @test mlc * Δ == MultiLeafCollimator(x.*Δ[1], y.*Δ[2])
    end
    
    y = -30.:5:30.

    n = length(y)-1
    x = rand(2, n)
    mlc = MultiLeafCollimator(x, y)
    
    @testset "Default Constructor" begin
        @test getpositions(mlc) == x
        @test getedges(mlc) == y
        @test length(mlc) == n
        @test size(mlc) == (n,)
    end

    @testset "Methods" begin

        @test firstindex(mlc) == 1
        @test lastindex(mlc) == n
        @test eachindex(mlc) == Base.OneTo(n)

        @test locate(mlc, -23.) == 2

        @test getedges(mlc) == y
        @test getedges(mlc, 2) == y[2:3]
        @test getedges(mlc, 2:5) == y[2:6]

        @test getpositions(mlc) == x
        @test getpositions(mlc, 3) == x[:, 3]
        @test getpositions(mlc, 3:6) == x[:, 3:6]

        xnew = rand(2, n)
        setpositions!(mlc, xnew)
        @test getpositions(mlc) == xnew

        # Similar
        mlc_sim = similar(mlc)
        @test length(mlc_sim) == length(mlc)
        @test getedges(mlc_sim) == getedges(mlc)

        # Deep copy
        mlc_cp = copy(mlc)
        @test mlc == mlc_cp
        setpositions!(mlc_cp, rand(size(mlc)))
        @test mlc != mlc_cp

        # Close Leaves
        closeleaves!(mlc_cp)
        @test all(getpositions(mlc_cp) .== 0.)
    end

    @testset "Indexing" begin
        i = 3
        @test mlc[i] == Jaws(mlc.positions[:, i], mlc.edges[i:i+1])

        i = 2:6
        x = mlc.positions[:, i]
        y = mlc.edges[i[1]:i[end]+1]
        @test mlc[i] == MultiLeafCollimator(x, y)

        i, xi = rand(1:n), rand(2)
        mlc[i] = xi
        @test mlc.positions[:, i] == xi

        i = 3:6
        xi = rand(2, length(i))
        mlc[i] = xi
        @test mlc.positions[:, i] == xi
        
    end

    @testset "Operate with Vector" test_operations(mlc, rand(2))
    @testset "Operate with Tuple" test_operations(mlc, tuple(rand(2)...))
        
    @testset "Basic Constructor" begin
        n = length(y)-1
        Δ = step(y)
        mlc = MultiLeafCollimator(n, Δ)
        
        @test getpositions(mlc) == zeros(2, n)

        x = rand(2, length(y)-1)
        setpositions!(mlc, x)
        @test mlc.positions == x
    end

    @testset "Zero Constructor" begin
        y = -60:10.:60
        mlc = MultiLeafCollimator(y)

        @test getedges(mlc) == y
        @test getpositions(mlc) == zeros(2, length(y)-1)
    end

    @testset "Assertion Errors" begin
        y = -30.:5:30.
        x = rand(2, length(y)-10)
        @test_throws AssertionError("Length of positions and edges do not match") MultiLeafCollimator(x, y)

        x = rand(1, length(y)-1)
        @test_throws AssertionError("Number of leaf positions per track != 2") MultiLeafCollimator(x, y)
    end

    @testset "REPL IO" begin
        @test length(Roentgen._str_closedaperture(10, 80)) == 80
        @test length(Roentgen._str_aperture(10, [10, 20], 80)) == 80
    end

end
