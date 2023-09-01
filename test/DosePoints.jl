@testset "Bounds" begin

    @testset "CylinderBounds" begin
        d, h = rand(2)
        c = SVector(rand(3)...)
        bounds = CylinderBounds(d, h, c)
    
        # Compare the two constructors
        @test CylinderBounds(d, h, zeros(3)) == CylinderBounds(d, h)
    
        pmin, pmax = Roentgen.extent(bounds)
        @test 0.5*(pmin+pmax) ≈ c
        @test pmax-pmin ≈ [d, d, h]
    
        # At center (inside)
        @test Roentgen.within(bounds, c)
    
        # Random position inside
        function rand_pos(α)
            θ = 2π*rand()
            @. c + 0.5*α*[d, d, h]*[cos(θ), sin(θ), 1.]
        end
        @test Roentgen.within(bounds, rand_pos(0.8))
    
        # Random position outside
        @test !Roentgen.within(bounds,  rand_pos(1.2))
    end

    @testset "SurfaceBounds" begin
        ϕ = 0.:π/4:2π
        y = -12.:2.:18.
        rhoval = 10.
        rho = fill(rhoval, length(ϕ), length(y))
        
        surf = CylindricalSurface(ϕ, y, rho, [0., 0., 0.])
        bounds = SurfaceBounds(surf)
        
        pmin, pmax = Roentgen.extent(bounds)
        @test pmin == [-rhoval, y[1], -rhoval]
        @test pmax == [ rhoval, y[end],  rhoval]

        # Inside
        
        @test Roentgen.within(bounds, zeros(3))
        
        function rand_pos(rho1, rho2, y1, y2)
            ϕi = 2π*rand()
            yi = (y2-y1)*rand()+y1
            rhoi = (rho2-rho1)*rand()+rho1
            [rhoi*cos(ϕi), yi, rhoi*sin(ϕi)]
        end
        
        p = rand_pos(0., rhoval, y[1], y[end])
        @test Roentgen.within(bounds, p)

        # Outside Radially
        p = rand_pos(rhoval, 2*rhoval, y[1], y[end])
        @test !Roentgen.within(bounds, p)
        
        # Outside Axially
        p = rand_pos(0, rhoval, y[1], y[1]-1)
        @test !Roentgen.within(bounds, p)

        p = rand_pos(0, rhoval, y[end], y[end]+1)
        @test !Roentgen.within(bounds, p)
    end
end

@testset "DosePoints" begin

    function test_operations(pos)
        @testset "$(String(Symbol(op)))" for op in (+, -, *, /)
            v = rand(3)
            pos_new = op(pos, v)
            
            for i=1:3
                @test getaxes(pos_new, i) == op.(getaxes(pos, i), v[i])
            end 
        end
    end

    @testset "DoseGrid" begin

        x, y, z = -10.:3.:10., 10.:5.:20, -20.:15.:30.
        axes = (x, y, z)
        pos = DoseGrid(axes)

        n = length.(axes)

        # Other Constructor
        @test pos == DoseGrid(x, y, z)

        # Base Methods
        
        @test size(pos) == n
        @test eachindex(pos) == CartesianIndices(n)
    
        @test Tuple(getaxes(pos)) == axes
        @test getaxes(pos, 2) == axes[2]
        
        # Indexing
        index = rand(CartesianIndices(n))
        @test pos[index] ≈ SVector(getindex.(axes, Tuple(index)))
        @test pos[Tuple(index)...] ≈ pos[index]

        test_operations(pos)

    end

    @testset "DoseGridMasked" begin

        x, y, z = -10.:3.:10., 10.:5.:20, -20.:15.:30.
        ax = x, y, z

        n = length.(ax)
        N = prod(n)
        nval = 5
        indices = getindex.(Ref(CartesianIndices(n)), sort(rand(1:N, nval)))
        
        pos = DoseGridMasked(ax, indices)

        @test size(pos) == (nval,)
        @test length(pos) == nval
        
        @test eachindex(pos) == Base.OneTo(nval)
        @test CartesianIndices(pos) == indices

        @test getaxes(pos) == ax
        @test getaxes(pos, 2) == ax[2]

        @testset "Indexing" for index in eachindex(pos)
            i, j, k = Tuple(indices[index])
            @test pos[index] == [x[i], y[j], z[k]]
        end
    
        # Operations
        test_operations(pos)

        # Constructor with bounds
        bounds = CylinderBounds(13., 12., SVector(rand(3)...))
        pos = DoseGridMasked(6., bounds)
        @test all(Roentgen.within.(Ref(bounds), pos))
        
    end
end
