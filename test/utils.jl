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

@testset "Scale to Isoplane" begin
    z_iso = rand()
    p = [2*rand().-1, 2*rand().-1, z_iso]
    
    @test DoseCalculations.scale_to_isoplane(p, z_iso) ≈ p[1:2]
    α = rand()
    @test DoseCalculations.scale_to_isoplane(p, α*z_iso) ≈ α*p[1:2]
    α = rand()+1
    @test DoseCalculations.scale_to_isoplane(p, α*z_iso) ≈ α*p[1:2]
end
