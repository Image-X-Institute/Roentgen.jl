
@testset "Dose Volume" begin
    x = -10:1.:10
    pos = DoseGrid(x, x, x)
    surf = PlaneSurface(800.)

    vol = DoseVolume(pos, surf)

    @test pos == getpositions(vol)
    @test surf == getsurface(vol)
end
