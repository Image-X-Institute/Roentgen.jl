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

end
