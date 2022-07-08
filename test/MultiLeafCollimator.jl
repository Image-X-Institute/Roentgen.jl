function test_subsetting(TVec)

    @testset "Subset Indices" begin
        y_test = TVec(-2.:1:2.)
        mlc = MultiLeafCollimator(y_test)

        @testset "Within leaf edges" begin
            y1, y2 = -0.8, 1.2
            iL, iU = DoseCalculations.subset_indices(mlc, y1, y2)
            @test iL == 2
            @test iU == 4
            @test y_test[iL]<=y1<y_test[iL+1]
            @test y_test[iU]<y2<=y_test[iU+1]
        end

        @testset "On lower leaf edge" begin
            y1, y2 = -1., 1.2
            iL, _ = DoseCalculations.subset_indices(mlc, y1, y2)
            @test iL == 2
            @test y_test[iL]<=y1<y_test[iL+1]
        end

        @testset "On upper leaf edge" begin
            y1, y2 = -0.8, 2.
            _, iU = DoseCalculations.subset_indices(mlc, y1, y2)
            @test iU == 4
            @test y_test[iU]<y2<=y_test[iU+1]
        end

        @testset "Below lower bound" begin
            y1, y2 = -2.8, 2.
            iL, _ = DoseCalculations.subset_indices(mlc, y1, y2)
            @test iL == 1
        end

        @testset "Above upper bound" begin
            y1, y2 = -0.8, 3.
            _, iU = DoseCalculations.subset_indices(mlc, y1, y2)
            @test iU == length(mlc)
        end
    end

    @testset "Subset of MLC" begin
        y_test = TVec(-10.:1:10.)
        mlc = MultiLeafCollimator(y_test)

        function test_subset(mlc, y1, y2)
            mlc_subset = extract_subset(mlc, y1, y2)
            @test edgeposition(mlc_subset, 1)<=y1
            @test y2<=edgeposition(mlc_subset, length(mlc_subset)+1)
        end

        @testset "Within leaf edges" begin
            test_subset(mlc, -4.5, 6.3)
        end

        @testset "On lower edge" begin
            test_subset(mlc, -5., 6.3) 
        end

        @testset "On upper edge" begin
            test_subset(mlc, -4.5, 6.) 
        end
        
    end
end

@testset "MultiLeafCollimator Tests" begin
    vector_types = [typeof(-1.:1.:1.), Vector{Float64}]
    @testset for TVec in vector_types
        test_subsetting(TVec)
    end
end