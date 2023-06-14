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

    @testset "SurfaceBounds" begin
        ϕ = deg2rad(-180):deg2rad(45.):deg2rad(180.)
        y = -12.:2.:18.
        ρval = 10.
        ρ = fill(ρval, length(ϕ), length(y))
        y1, y2 = y[[1,end]]
        
        surf = CylindricalSurface(ϕ, y, ρ)
        bounds = SurfaceBounds(surf)
        
        pmin, pmax = DoseCalculations.extent(bounds)
        @test pmin == [-ρval, y1, -ρval]
        @test pmax == [ ρval, y2,  ρval]

        # Inside
        
        @test DoseCalculations.within(bounds, zeros(3))
        
        function rand_pos(ρlim, ylim)    
            ϕ = 2π*rand()
            y = (ylim[end]-ylim[1])*rand()+ylim[1]
            ρ = (ρlim[end]-ρlim[1])*rand()+ρlim[1]
            [ρ*cos(ϕ), y, ρ*sin(ϕ)]
        end
        
        @test DoseCalculations.within(bounds, rand_pos([0, 1]*ρval, y))

        # Outside Radially
        @test !DoseCalculations.within(bounds, rand_pos([1, 2]*ρval, y))
        
        # Outside Axially
        ϕi, ρi = (2π, ρval).*rand(2)
        x, z = ρi*cos(ϕi), ρi*sin(ϕi)
        @test !DoseCalculations.within(bounds, [x, y1-1., z])

        ϕi, ρi = (2π, ρval).*rand(2)
        x, z = ρi*cos(ϕi), ρi*sin(ϕi)
        @test !DoseCalculations.within(bounds, [x, y2+1., z])
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
        @test all(DoseCalculations.within.(Ref(bounds), pos))
        
    end
end
