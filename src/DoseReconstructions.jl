"""
    Dose Reconstructions

High-level functions to compute dose reconstructions from a given set of beams.
"""

export reconstruct_dose

"""
    reconstruct_dose(pos, surf, plan, calc; Δx=1., show_progess=true)

Reconstruct the dose.

Requires dose positions (`pos`), an external surface (`surf`), a set of beams
(`plan`), and the dose calculation algorithm (`calc`). Optional arguments
include:
- `Δx`: Size of each bixel in the fluence grid (defaults to 5.)
- `ΔMU`: Meterset increments (defaults to 2.)
- `show_progess`: Whether to display the progress (defaults to `true`)
"""
function reconstruct_dose(pos, surf, plan, calc; Δx=5., show_progress=true)
    dose = zeros(length(pos))
    Ψ = zeros(1)

    # Iterate through each field in the plan
    for field in plan

        pos_fixed = patient_to_fixed(getisocenter(field)).(pos)

        if(show_progress)
            totalMU = 0.
            p = Progress(length(field))
        end

        # Iterate through each control point in the field
        for beam in field
            ΔMU = getΔMU(beam)
            gantry = getgantry(beam)
            mlc = getmlc(beam)
            jaws = getjaws(beam)

            bixels = bixels_from_bld(mlc, jaws; Δx=Δx)
            beamlets = Beamlet.(bixels, Ref(gantry))
            nbixels = length(bixels)

            resize!(Ψ, nbixels)
        
            dose .+= ΔMU*sum(dose_fluence_matrix(SparseMatrixCSC, vec(pos_fixed),
                                                 vec(beamlets), surf, calc);
                             dims=2)
            
            if(show_progress)
                next!(p; showvalues=[(:MU, totalMU+=ΔMU)])
            end

        end
    end
    dose
end
