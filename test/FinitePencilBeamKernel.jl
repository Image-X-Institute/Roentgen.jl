@testset "Pencil Beam Profile" begin
    # Tests if profile is continuous
    @testset "Continuity" begin
        u, x₀ = rand(2)
        @test DoseCalculations.fpkb_profile_left(-x₀, u, x₀) ≈ DoseCalculations.fpkb_profile_center(-x₀, u, x₀)
        @test DoseCalculations.fpkb_profile_right(x₀, u, x₀) ≈ DoseCalculations.fpkb_profile_center(x₀, u, x₀)
    end

    # Tests Eq. 4 in paper
    @testset "1D Profile Integration" begin
        u, x₀ = rand(2)
        F = x->DoseCalculations.fpkb_profile(x, u, x₀)
        I, ε = quadgk(F, -Inf, Inf) # quadgk returns integral and upper bound of error
        @test norm(I - 2*x₀) < ε
    end

    # Tests Eq. 7 in paper
    @testset "2D Profile Integration" begin
        
        w = rand()
        ux = [rand(), rand()]
        uy = [rand(), rand()]
        x₀ = rand()
        y₀ = rand()

        F = (x)->DoseCalculations.fpkb_dose(x[1], x[2], w, ux, uy, x₀, y₀)
        
        lb, ub = -1e3, 1e3
        I, ε = hcubature(F, fill(lb, 2), fill(ub, 2))
        @test norm(I - 4*x₀*y₀) < ε
    end

end
