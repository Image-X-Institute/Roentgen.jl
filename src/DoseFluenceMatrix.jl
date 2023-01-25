export dose_fluence_matrix, dose_fluence_matrix!

#--- Dose-Fluence Matrix-------------------------------------------------------

"""
    dose_fluence_matrix(pos, bixels, bixels, surf, calc)

Compute a dose-fluence matrix from dose positions, bixels, external surface and
dose calculation algorithm.

See `dose_fluence_matrix!` for implementation.
"""
function dose_fluence_matrix(pos, bixels::AbstractVector{<:AbstractFluenceElement},
                             gantry::GantryPosition,
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)
    D = spzeros(length(bixels), length(pos))
    dose_fluence_matrix!(D, pos, bixels, gantry, surf, calc)
end

"""
    dose_fluence_matrix!(D, pos, bixels, bixels, surf, calc)

Compute a dose-fluence matrix from dose positions, bixels, external surface and
dose calculation algorithm.

Requires the `point_kernel!` method to be defined for the given dose calculation
algorithm (`calc`). `point_kernel!` computes the dose calculated from the set of
bixels a given dose point. Stores result in `D`.
"""
function dose_fluence_matrix!(D::SparseMatrixCSC, pos, bixels::AbstractVector{<:AbstractFluenceElement},
                              gantry::GantryPosition,
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)

    colptr = D.colptr
    rowval = D.rowval
    nzval = D.nzval

    # Fill colptr
    colptr[1] = 1
    @batch per=thread for j in eachindex(pos) #
         # Temporarily store number of values in colptr
        colptr[j+1] = kernel_size(calc, pos[j], bixels, gantry)
    end
    cumsum!(colptr, colptr) # Compute colptr from number of values

    # Preallocate arrays
    n_prealloc = colptr[end]-1
    resize!(nzval, n_prealloc)
    resize!(rowval, n_prealloc)

    # # Compute row and matrix values
    @batch per=thread for j in eachindex(pos) #

        ptr = colptr[j]:(colptr[j+1]-1)

        # Create views of rowval and nzval
        I = @view rowval[ptr]
        V = @view nzval[ptr]
        
        # Calculate the row and matrix value for given position
        point_kernel!(I, V, pos[j], bixels, surf, gantry, calc)
    end

    dropzeros!(D)

    D
end
