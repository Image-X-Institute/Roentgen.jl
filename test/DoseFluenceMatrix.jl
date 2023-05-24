

@testset "Kernel Truncation" begin
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
