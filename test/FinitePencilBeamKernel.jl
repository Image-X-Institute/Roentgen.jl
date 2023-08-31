@testset "Pencil Beam Profile" begin
    # Tests if profile is continuous
    @testset "Continuity" begin
        u, x₀ = [0.5, 7.]
        @test Roentgen.fpbk_profile_left(-x₀, u, x₀) ≈ Roentgen.fpbk_profile_center(-x₀, u, x₀)
        @test Roentgen.fpbk_profile_right(x₀, u, x₀) ≈ Roentgen.fpbk_profile_center(x₀, u, x₀)
    end

    # Tests Eq. 4 in paper
    @testset "1D Profile Integration" begin
        u, x₀ = [0.5, 7.]
        F = x->Roentgen.fpbk_profile(x, u, x₀)
        I, ε = quadgk(F, -Inf, Inf) # quadgk returns integral and upper bound of error
        @test norm(I - 2*x₀) < ε
    end

    # Tests Eq. 7 in paper
    @testset "2D Profile Integration" begin
        
        w = 0.3
        ux = [0.5, 0.04]
        uy = [0.4, 0.03]
        x₀, y₀ = [7., 6.]
        F = (x)->Roentgen.fpbk_dose(x[1], x[2], w, ux, uy, x₀, y₀)
        
        lb, ub = -1e3, 1e3
        integral, ε = hcubature(F, fill(lb, 2), fill(ub, 2))

        @test norm(integral - 4*x₀*y₀) < ε
    end

end

@testset "Finite Pencil Beam Kernel" begin

    linear(x, a, b) = a*x+b
    
    function setup_calc(wfun, uxfun, uyfun, Afun)
        depth = range(0., 300., length=11)
    
        w = wfun.(depth)
        ux = uxfun.(depth)
        uy = uyfun.(depth)
    
        parameters = [SVector(w[i], ux[i]..., uy[i]...) for i in eachindex(w, ux, uy)]
    
        tanθ = range(0., deg2rad(1.), length=9)
    
        scalingfactor = @. Afun(depth, tanθ')
    
        depth, tanθ, FinitePencilBeamKernel(parameters, scalingfactor, depth, tanθ)
    end

    wfun(x) = linear(x, 0.001, 0.4)
    uxfun(x) = [0.5, 0.04].*linear(x, 0.001, 0.4)
    uyfun(x) = [0.4, 0.03].*linear(x, 0.001, 0.4)
    Afun(x, y) = linear(x, 0.001, 0.4)*linear(y, 0.1, 0.1)
    
    depth, tanθ, calc = setup_calc(wfun, uxfun, uyfun, Afun)

    @testset "Parameters and Scaling Factor" begin
        
        d = depth[end]*rand()
        w, ux, uy = Roentgen.getparams(calc, d) 
                
        @test w ≈ wfun(d)
        @test ux ≈ uxfun(d)
        @test uy ≈ uyfun(d)
        
        d = depth[end]*rand()
        t = tanθ[end]*rand()
        
        A = Roentgen.getscalingfactor(calc, d, t)
        @test A ≈ Afun(d, t)
    end

    @testset "Depth Check" begin
        T = Float64
        @test !Roentgen._depth_check(zero(T))
        @test !Roentgen._depth_check(eps(T))
        @test Roentgen._depth_check(-eps(T))
        @test Roentgen._depth_check(Inf)
        @test Roentgen._depth_check(NaN)
    
        λ = 100.
        surf = ConstantSurface(λ)
        gantry = GantryPosition(2π*rand(), 0., λ)
        beamlet = Beamlet(Bixel(0., 0., 1., 1.), gantry)

        src = Roentgen.getposition(gantry)

        posfun(λ) = src - λ*normalize(src)
        @test Roentgen.point_dose(posfun(1.1*λ), beamlet, surf, calc) > 0.
        @test Roentgen.point_dose(posfun(λ), beamlet, surf, calc) > 0.
        @test Roentgen.point_dose(posfun(0.9*λ), beamlet, surf, calc) ≈ 0.

    end

end
