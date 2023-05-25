

@testset "Kernel Truncation" begin

    # Pre-calculated tests
    a = SVector(0., 0., 1.)
    SAD = 2.
    maxradius = 3.
    @testset "Pre-calculated Tests" begin
        r_inside = [SVector(1.,0.,1.), SVector(-2.,0.,3.), SVector(3.,0.,4.)]
        r_outside = [SVector(4.,0.,2.), SVector(2.,0.,1.), SVector(-6.,0.,3.)]
        for r in r_inside
            @test DoseCalculations.kernel_size(r, a, maxradius/SAD)
        end
        for r in r_outside
            @test !DoseCalculations.kernel_size(r, a, maxradius/SAD)
        end

        for r in r_inside
            trans = RotYZX(2*π*rand(3)...)
            @test DoseCalculations.kernel_size(trans*r, trans*a, maxradius/SAD)
        end
        for r in r_outside
            trans = RotYZX(2*π*rand(3)...)
            @test !DoseCalculations.kernel_size(trans*r, trans*a, maxradius/SAD)
        end
    end

    # Arbitrary a

    rand_centred(args...) = 2*rand(args...).-1

    random_vector() = normalize(SVector(rand_centred(3)...))

    # Any r along a will return true
    a = random_vector()
    r = 2*rand()*a
    @test DoseCalculations.kernel_size(r, a, rand())

    random_perp_vector(a) = normalize(cross(a, rand_centred(3)))
    function rand_setup(α)
        # Create a random beamlet direction
        a = random_vector()
        # Create a random vector perpendicular to a
        δ = random_perp_vector(a)

        # Random maxradius and SAD
        maxradius = 2*rand()
        SAD = 2*rand()

        # Set r at a perpendicular distance away from a as isoplane, then scale
        # by a random amount
        r = 2*rand()*normalize(SAD*a + α*maxradius*δ)

        a, r, SAD, maxradius
    end

    # Test inside radius
    a, r, SAD, maxradius = rand_setup(0.9)
    @test DoseCalculations.kernel_size(r, a, maxradius/SAD)

    # Test outside radius
    a, r, SAD, maxradius = rand_setup(1.1)
    @test !DoseCalculations.kernel_size(r, a, maxradius/SAD)

end
