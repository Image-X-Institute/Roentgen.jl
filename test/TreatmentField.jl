
function rand_mlc(y)
    n = length(y)-1
    xB = rand(n)
    xA = xB .+ rand(n)
    x = vcat(xB', xA')
    @. x = 20*x-10.
    
    MultiLeafCollimator(x, y)
end

@testset "Control Point" begin
    mlc = rand_mlc(-15.:5.:15)
    jaws = Jaws(-13., 14., -23., 22.)
    ϕg, θb, SAD, ΔMU, dose_rate = rand(5)
    isocenter = SVector(rand(3)...)

    pt = DoseCalculations.ControlPoint(mlc, jaws, ϕg, θb, SAD, ΔMU, dose_rate, isocenter)

    @testset "Methods" begin
        @test getmlc(pt) == mlc
        @test getjaws(pt) == jaws
        @test getdoserate(pt) == dose_rate
        @test getisocenter(pt) == isocenter

        @test fixed_to_bld(pt) == fixed_to_bld(ϕg, θb, SAD)
        @test getgantry(pt) == GantryPosition(ϕg, θb, SAD)
        @test getΔMU(pt) == ΔMU
    end
end

@testset "VMATField" begin

    ncontrol = 5
    y = -15:5.:15
    x = rand(2, length(y)-1, ncontrol)
    mlc = MultiLeafCollimatorSequence(x, y);
    jaws = Jaws(-13., 14., -23., 22.)
    
    ϕg = rand(ncontrol)
    meterset = cumsum(rand(ncontrol))
    meterset .-= meterset[1]
    θb, SAD, dose_rate = rand(3)
    isocenter = SVector(rand(3)...)
    
    field = VMATField(mlc, jaws, ϕg, θb, SAD, meterset, dose_rate, isocenter)
    
    @test length(field) == ncontrol
    
    @testset "Methods" begin
        @test getmlc(field) == mlc
        @test getjaws(field) == jaws
        @test getdoserate(field) == dose_rate
        @test getisocenter(field) == isocenter
    
        @test getgantry(field) == GantryPosition.(ϕg, θb, SAD)
        @test getmeterset(field) == meterset
    end
    
    @testset "Index Methods" begin
        i = rand(1:ncontrol)
        @test getmlc(field, i) == mlc[i]
        @test getjaws(field, i) == jaws
        @test getdoserate(field, i) == dose_rate
        @test getisocenter(field, i) == isocenter
    
        @test getgantry(field, i) == GantryPosition(ϕg[i], θb, SAD)
        @test getmeterset(field, i) == meterset[i]
    
        @test getΔMU(field, 1) == 0.5*(meterset[2]-meterset[1])
        @test getΔMU(field, 3) == 0.5*(meterset[4]-meterset[2])
        @test getΔMU(field, ncontrol) == 0.5*(meterset[ncontrol]-meterset[ncontrol-1])
    
        @test sum(getΔMU.(Ref(field), 1:ncontrol)) ≈ meterset[end]
    end

    @testset "Resample" begin
        field = VMATField(mlc, jaws, ϕg, θb, SAD, meterset, dose_rate, isocenter)
        field_r = resample(field, 0.1)

        @testset "$(Symbol(f))" for f in [getgantry, getmlc, getjaws, getmeterset]
            @test f(field_r, 1) == f(field, 1)
            @test f(field_r, length(field_r)) == f(field, ncontrol)
        end
    end

end

@testset "Fix Angle" begin
    @test DoseCalculations.fixangle(0.4π) == 0.4π
    @test DoseCalculations.fixangle(π) == π
    @test DoseCalculations.fixangle(1.3π) == (1.3-2)π
end
