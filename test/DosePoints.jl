@testset "DosePoints" begin
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

        function test_hdf5(pos)
            filename = "tmp.hdf5"
            
            h5open(filename, "w") do file
                save(file, pos)
            end

            @testset "Write" begin 
                h5open(filename, "r") do file
                    @test haskey(file, "pos")
                    @test haskey(file, "pos/x")
                    @test haskey(file, "pos/y")
                    @test haskey(file, "pos/z")
                end
            end

            @testset "Read" begin 
                pos2 = h5open(filename, "r") do file
                    load(DoseGrid, file)
                end
                @test pos2.axes == pos.axes
            end

            rm(filename)
        end

        pos = DoseGrid(5., CylinderBounds(200., 200.))
        test_operations(pos)
        @testset "IO - HDF5" test_hdf5(pos)

    end

    @testset "DoseGridMasked" begin
        
        function test_hdf5(pos)
            filename = "tmp.hdf5"
            
            h5open(filename, "w") do file
                save(file, pos)
            end

            @testset "Write" begin 
                h5open(filename, "r") do file
                    @test haskey(file, "pos")
                    @test haskey(file, "pos/x")
                    @test haskey(file, "pos/y")
                    @test haskey(file, "pos/z")
                    @test haskey(file, "cells/index")
                    @test haskey(file, "cells/cells")
                end
            end

            @testset "Read" begin 
                pos2 = h5open(filename, "r") do file
                    load(DoseGridMasked, file)
                end
                @test pos2.axes == pos.axes
            end
            
            rm(filename)
        end

        pos = DoseGridMasked(5., CylinderBounds(200., 200.))
        test_operations(pos)
        @testset "IO - HDF5" test_hdf5(pos)
    end
end