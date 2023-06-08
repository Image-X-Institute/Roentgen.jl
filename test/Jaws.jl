@testset "Jaws" begin

    x = rand(2)
    y = rand(2)
    jaws = Jaws(x, y)
    
    @test getx(jaws) == x
    @test gety(jaws) == y
    
    @test jaws == Jaws(x[1], x[2], y[1], y[2])
    
    fieldsize = rand()
    jaws = Jaws(fieldsize)
    @test getx(jaws) == 0.5*fieldsize*[-1, 1]
    @test gety(jaws) == 0.5*fieldsize*[-1, 1]

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
