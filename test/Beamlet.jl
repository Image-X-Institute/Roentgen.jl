@testset "Beamlet" begin

    SAD = 100*rand()

    pos = 20*rand(2).-1
    width = 2*rand(2)

    bixel = Bixel(pos, width)
    gantry = GantryPosition(2π*rand(), π*rand(), SAD)

    beamlet = Beamlet(bixel, gantry)

    # Source Axis
    s = Roentgen.beamaxis(gantry)
    @test Roentgen.source_axis(beamlet) ≈ s
    @test Roentgen.source_axis_distance(beamlet) == SAD
    @test Roentgen.source_position(beamlet) == SAD*s

    # Beamlet Axes
    axes = Roentgen.beamlet_axes(beamlet)

    # - All be unit vectors
    @test all(@. norm(axes) ≈ 1)
    # - All be perpendicular
    tol = 1e-12
    @test isapprox(dot(axes[1], axes[2]), 0., atol=tol)
    @test isapprox(dot(axes[2], axes[3]), 0., atol=tol)
    @test isapprox(dot(axes[3], axes[1]), 0., atol=tol)

    @test Roentgen.direction(beamlet) == axes[3]

    @test Roentgen.halfwidth(beamlet) ≈ 0.5*width
    @test Roentgen.width(beamlet) ≈ width

end
