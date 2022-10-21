@testset "VMATField" begin
    
    @testset "HDF5" begin

        n = 10

        y = -30:5.:30
        x = rand(2, length(y)-1, n)
        mlc = MultiLeafCollimatorSequence(x, y)

        jaws = Jaws(rand(2), rand(2))

        ϕg = rand(n)
        θb = rand()
        SAD = rand()

        meterset = rand(n)
        Ḋ = rand()
        isocenter = rand(3)

        field = VMATField(mlc, jaws, ϕg, θb, SAD, meterset, Ḋ, isocenter)

        filename = "tmp.hdf5"
        h5open(filename, "w") do file
            save(file, field)
        end

        @testset "Write/Read" begin
            field2 = h5open(filename, "r") do file
                load(VMATField, file)
            end
            @test field.mlc == field2.mlc
            @test field.jaws == field2.jaws
            @test field.gantry_angle == field2.gantry_angle
            @test field.collimator_angle == field2.collimator_angle
            @test field.source_axis_distance == field2.source_axis_distance
            @test field.meterset == field2.meterset
            @test field.dose_rate == field2.dose_rate
            @test field.isocenter == field2.isocenter
        end

        rm(filename)

    end
end
