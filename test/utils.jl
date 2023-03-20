@testset "Snapped Range" begin

    Δ = 1.
    x1 = -10.5
    x2 = 10.5

    function test_snapped_range(x1, x2, Δ)
        r = DoseCalculations.snapped_range(x1, x2, Δ)
        @test r[1]<=x1
        @test x2<=r[end]
    end

    test_snapped_range(-10.8, 10.8, 1.)
    test_snapped_range(-10., 10., 1.)
    test_snapped_range(-10*rand(), 10*rand(), 2*rand())

end