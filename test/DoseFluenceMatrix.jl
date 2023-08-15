

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

@testset "Dose-Fluence Matrix Calculation" begin
    depth = 0.:20.:400
    x = -200.:20.:200.
    pos = SVector.(x, 0., -depth')

    fs = 10.
    xb = -0.5*fs:5:0.5*fs
    bixels = BixelGrid(xb, xb)

    calc = MockKernel()

    surf = PlaneSurface(1000.)

    gantry = GantryPosition(0., 0., 1000.)
    beamlets = Beamlet.(bixels, (gantry,))

    Dsaved = JLD2.load("test-data/dose_fluence_matrix.jld2", "mock-kernel")

    @test Dsaved == dose_fluence_matrix(Matrix, vec(pos), vec(beamlets), surf, calc; maxradius=5.)
    @test Dsaved == dose_fluence_matrix(SparseMatrixCSC, vec(pos), vec(beamlets), surf, calc; maxradius=5.)
end
