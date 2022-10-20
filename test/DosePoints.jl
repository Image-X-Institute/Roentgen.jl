
function test_operations(pos::DoseCalculations.AbstractDoseGrid)
    @testset "Operations" begin
        @testset "$(String(Symbol(op)))" for op in (+, -, *, /)
            v = rand(3)
            pos_new = op(pos, v)
            
            @test getx(pos_new) == op.(getx(pos), v[1])
            @test gety(pos_new) == op.(gety(pos), v[2])
            @test getz(pos_new) == op.(getz(pos), v[3])    
        end
    end
end

@testset "DoseGrid" begin
    pos = DoseGrid(5., CylinderBounds(200., 200.))
    test_operations(pos)
end

@testset "DoseGridMasked" begin
    pos = DoseGridMasked(5., CylinderBounds(200., 200.))
    test_operations(pos)
end