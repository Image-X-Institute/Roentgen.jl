@testset "Pencil Beam Profile" begin
    # Tests if profile is continuous
    @testset "Continuity" begin
        u, x₀ = [0.5, 7.]
        @test DoseCalculations.fpbk_profile_left(-x₀, u, x₀) ≈ DoseCalculations.fpbk_profile_center(-x₀, u, x₀)
        @test DoseCalculations.fpbk_profile_right(x₀, u, x₀) ≈ DoseCalculations.fpbk_profile_center(x₀, u, x₀)
    end

    # Tests Eq. 4 in paper
    @testset "1D Profile Integration" begin
        u, x₀ = [0.5, 7.]
        F = x->DoseCalculations.fpbk_profile(x, u, x₀)
        I, ε = quadgk(F, -Inf, Inf) # quadgk returns integral and upper bound of error
        @test norm(I - 2*x₀) < ε
    end

    # Tests Eq. 7 in paper
    @testset "2D Profile Integration" begin
        
        w = 0.3
        ux = [0.5, 0.04]
        uy = [0.4, 0.03]
        x₀, y₀ = [7., 6.]
        F = (x)->DoseCalculations.fpbk_dose(x[1], x[2], w, ux, uy, x₀, y₀)
        
        lb, ub = -1e3, 1e3
        integral, ε = hcubature(F, fill(lb, 2), fill(ub, 2))

        @test norm(integral - 4*x₀*y₀) < ε
    end

end
