rand_in_range(x1, x2) = (x2-x1)*rand() + x1

function test_grid(xg)
    @testset "Inbounds" begin
        xi = rand_in_range(xg[1], xg[end])
        i = locate(xg, xi)
        @test xg[i]<=xi<xg[i+1]
    end

    @testset "Out of bounds xi<xg[1]" begin
        xi = xg[1] - rand()
        i = locate(xg, xi)
        @test i<=1
    end

    @testset "Out of bounds xi>=xg[end]" begin
        xi = xg[end] + rand()
        i = locate(xg, xi)
        @test length(xg)-1<=i
    end
end

@testset "Location" begin
    @testset "Uniform Spacing" begin
        xg = range(0., 1., length=11)
        test_grid(xg)
    end
    @testset "Non-Uniform Spacing" begin
        n = 10
        xg = range(0., 1., length=n+1) .+ rand_in_range(-1., 1.)*0.5/n
        test_grid(xg)
    end
end

function test_linear_interpolation(xg)
    fg = rand(2)

    xs = [0., 1.]

    @testset "f[$(i[1])] == f($(xs[i[1]]))" for i in eachindex(fg)
        @test fg[i] == interp(xg, fg, xs[i])
    end

    @testset "avg(f) == f(0.5,0.5)" begin
        @test sum(fg)/length(fg) == interp(xg, fg, 0.5)
    end
end

@testset "Linear Interpolation" begin
    @testset "Uniform Grid" begin
        xg = 0.:1.:1.
        test_linear_interpolation(xg)
    end
    @testset "Non-uniform Grid" begin
        xg = [0., 1.]
        test_linear_interpolation(xg)
    end
end

function test_bilinear_interpolation(xg)
    yg = copy(xg)
    fg = rand(2, 2)

    xs = [0., 1.]
    ys = [0., 1.]

    @testset "f[$(i[1]),$(i[2])] == f($(xs[i[1]]),$(ys[i[2]]))" for i in CartesianIndices(fg)
        xi = xs[i[1]]
        yi = ys[i[2]]
        fi = interp(xg, yg, fg, xi, yi)
        @test fi == fg[i]
    end

    @testset "avg(f) == f(0.5,0.5)" begin
        @test sum(fg)/length(fg) == interp(xg, yg, fg, 0.5, 0.5)
    end
end

@testset "Bilinear Interpolation" begin
    @testset "Uniform Grid" begin
        xg = 0.:1.:1.
        test_bilinear_interpolation(xg)
    end
    @testset "Non-uniform Grid" begin
        xg = [0., 1.]
        test_bilinear_interpolation(xg)
    end
end

