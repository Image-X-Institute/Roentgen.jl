using HDF5

@testset "HDF5" begin

    jaws = Jaws(rand(2), rand(2))

    filename = "tmp.hdf5"
    h5open(filename, "w") do file
        save(file, jaws)
    end

    @testset "Write" begin
        h5open(filename, "r") do file
            @test haskey(file, "jaws")
            @test haskey(file, "jaws/x")
            @test haskey(file, "jaws/y")

            @test read(file["jaws/x"]) == Vector(getx(jaws))
            @test read(file["jaws/y"]) == Vector(gety(jaws))
        end

    end

    @testset "Read" begin
        jaws2 = h5open(filename, "r") do file
            load(Jaws, file)
        end
        @test getx(jaws2) == getx(jaws)
        @test gety(jaws2) == gety(jaws)
    end

    rm(filename)

end
