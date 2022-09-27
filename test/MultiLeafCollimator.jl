begin "MultiLeafCollimator"

    function test_constructor(mlc, x, y)
        @test mlc.positions == x
        @test mlc.edges == y
        @test mlc.n == length(y)-1
        @test mlc.n == size(x, 2)
    end

    function test_constructor_errors(mlc, x, y)

        @test mlc.positions == x
        @test mlc.edges == y
        @test mlc.n == length(y)-1
        @test mlc.n == size(x, 2)
    end

    function test_methods(mlc)
        @test length(mlc) == mlc.n
        @test size(mlc) == (mlc.n,)
        @test firstindex(mlc) == 1
        @test lastindex(mlc) == mlc.n
        @test eachindex(mlc) == Base.OneTo(mlc.n)
    end

    function test_indexing(mlc)
        i = 3
        @test mlc[i] == (mlc.positions[:, i], mlc.edges[i:i+1])

        i = 2:6
        x = mlc.positions[:, i]
        y = mlc.edges[i[1]:i[end]+1]
        @test mlc[i] == MultiLeafCollimator(x, y)
    end

    function test_operations(mlc, Δ)
        x = mlc.positions
        y = mlc.edges
        
        @test mlc + Δ == MultiLeafCollimator(x.+Δ[1], y.+Δ[2])
        @test mlc - Δ == MultiLeafCollimator(x.-Δ[1], y.-Δ[2])
        @test mlc / Δ == MultiLeafCollimator(x./Δ[1], y./Δ[2])
        @test mlc * Δ == MultiLeafCollimator(x.*Δ[1], y.*Δ[2])
    end

    y = -30.:5:30.
    @testset "With $(typeof(yᵢ)) edges" for yᵢ in [y, collect(y)]
        x = rand(2, length(yᵢ)-1)
        mlc = MultiLeafCollimator(x, yᵢ)
        
        @testset "Constructor" test_constructor(mlc, x, yᵢ)

        @testset "Constructor Errors" test_constructor(mlc, x, yᵢ)

        @testset "Basic Methods" test_methods(mlc)
        @testset "Indexing" test_indexing(mlc)
        
        @testset "Operate with Vector" test_operations(mlc, rand(2))
        @testset "Operate with Tuple" test_operations(mlc, tuple(rand(2)...))
    end

    @testset "Basic Constructor" begin
        y = -60:10.:60
        n = length(y)-1
        Δ = step(y)
        mlc = MultiLeafCollimator(n, Δ)
        
        @test mlc.edges == y
        @test mlc.positions == zeros(2, n)
    end

    y = -30.:5:30.
    @testset "Assertion Errors" begin
        x = rand(2, length(y)-10)
        @test_throws AssertionError("Length of positions and edges do not match") MultiLeafCollimator(x, y)

        x = rand(1, length(y)-1)
        @test_throws AssertionError("Number of leaf positions per track != 2") MultiLeafCollimator(x, y)
    end
end