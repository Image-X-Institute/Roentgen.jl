
@testset "Dose Volume" begin
    x = -10:10.:10
    pos = DoseGrid(x, x, x)
    surf = ConstantSurface(800.)

    vol = DoseVolume(pos, surf)

    @test pos == getpositions(vol)
    @test surf == getsurface(vol)

    beamlet = Beamlet(Bixel(0., 0., 1., 1.), GantryPosition(0., 0., 1000.))
    calc = MockKernel()

    D_truth = dose_fluence_matrix(Matrix, vec(pos), [beamlet], surf, calc)

    D_vol = dose_fluence_matrix(Matrix, vol, [beamlet], calc)
    @test D_vol ≈ D_truth

    D_vol .= 0
    dose_fluence_matrix!(D_vol, vol, [beamlet], calc)
    @test D_vol ≈ D_truth
    
end
