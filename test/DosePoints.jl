@testset "Bounds" begin

    @testset "CylinderBounds" begin
        d, h = rand(2)
        c = SVector(rand(3)...)
        bounds = CylinderBounds(d, h, c)
    
        # Compare the two constructors
        @test CylinderBounds(d, h, zeros(3)) == CylinderBounds(d, h)
    
        pmin, pmax = DoseCalculations.extent(bounds)
        @test 0.5*(pmin+pmax) ≈ c
        @test pmax-pmin ≈ [d, d, h]
    
        # At center (inside)
        @test DoseCalculations.within(bounds, c)
    
        # Random position inside
        function rand_pos(α)
            θ = 2π*rand()
            @. c + 0.5*α*[d, d, h]*[cos(θ), sin(θ), 1.]
        end
        @test DoseCalculations.within(bounds, rand_pos(0.8))
    
        # Random position outside
        @test !DoseCalculations.within(bounds,  rand_pos(1.2))
    end
end

@testset "DosePoints" begin

    function test_operations(pos::DoseCalculations.AbstractDoseGrid)
        @testset "Operations" begin
            @testset "$(String(Symbol(op)))" for op in (+, -, *, /)
                v = rand(3)
                pos_new = op(pos, v)
                
                for i=1:3
                    @test getaxes(pos_new, i) == op.(getaxes(pos, i), v[i])
                end 
            end
        end
    end

    @testset "DoseGrid" begin

        axes = (-10.:1.:10., 10.:2.:20, -20.:5.:30.)
        pos = DoseGrid(axes...)
        
        n = length.(axes)
        N = prod(n)
        
        @test size(pos) == n
        @test length(pos) == N
        
        @test eachindex(pos) == Base.OneTo(N)
        @test CartesianIndices(pos) == CartesianIndices(n)

        @test Tuple(getaxes(pos)) == axes
        
        # Linear Indexing
        index = rand(1:N)
        gridindex = Tuple(CartesianIndices(n)[index])
        @test pos[index] ≈ SVector(getindex.(axes, gridindex))
        
        # Cartesian Indexing
        index = rand(CartesianIndices(n))
        @test pos[index] ≈ SVector(getindex.(axes, Tuple(index)))
        @test pos[Tuple(index)...] ≈ pos[index]
        
        # Iteration
        p, i = iterate(pos)
        @test p ≈ pos[1] && i==2
        
        p, i = iterate(pos, 2)
        @test p ≈ pos[2] && i==3

        test_operations(pos)

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

        @testset "IO - HDF5" test_hdf5(pos)

    end

    @testset "DoseGridMasked" begin

        axes = (-10.:1.:10., 10.:2.:20, -20.:5.:30.)

        n = length.(axes)
        N = prod(n)
        
        nval = 5
        indices = getindex.(Ref(CartesianIndices(n)), sort(rand(1:N, nval)))
        
        pos = DoseGridMasked(axes..., indices, Vector{Vector{Int}}(undef, 0))
        
        @test size(pos) == (nval,)
        @test length(pos) == nval
        
        @test eachindex(pos) == Base.OneTo(nval)
        @test CartesianIndices(pos) == indices

        @test Tuple(getaxes(pos)) == axes
        
        # Linear Indexing
        index = rand(1:nval)
        gridindex = CartesianIndices(pos)[index]
        @test pos[index] ≈ SVector(getindex.(axes, Tuple(gridindex)))
        
        # Cartesian Indexing
        index = rand(CartesianIndices(pos))
        @test pos[index] ≈ SVector(getindex.(axes, Tuple(index)))
        @test pos[Tuple(index)...] ≈ pos[index]
        
        # Iteration
        p, i = iterate(pos)
        @test p ≈ pos[1] && i==2
        
        p, i = iterate(pos, 2)
        @test p ≈ pos[2] && i==3

        # Operations
        test_operations(pos)

        # Constructor with bounds
        bounds = CylinderBounds(13., 12., SVector(rand(3)...))
        pos = DoseGridMasked(2., bounds)
        @test all([DoseCalculations.within(bounds, pos[i]) for p in pos])
        
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

        @testset "IO - HDF5" test_hdf5(pos)
    end
end
