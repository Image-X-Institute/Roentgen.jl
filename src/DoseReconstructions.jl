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
- `Δx`: Size of each bixel in the fluence grid (defaults to 1.)
- `ΔMU`: Meterset increments (defaults to 2.)
- `show_progess`: Whether to display the progress (defaults to `true`)
"""
function reconstruct_dose(pos, surf, plan, calc; Δx=1., show_progress=true)
    # Allocate dose arrays
    dose = zeros(length(pos))


    # Iterate through each field in the plan
    for field in plan

        beams = field #resample(field, ΔMU, by=:MU)
        
        bixels = vec(bixel_grid(getmlc(field, 1), getjaws(field), Δx))
        D = spzeros(length(pos), length(bixels))
        Ψ = zeros(size(bixels))

        pos_fixed = patient_to_fixed(getisocenter(field)).(pos)

        if(show_progress)
            totalMU = 0.
            p = Progress(length(beams))
        end

        beamlets = Vector{Beamlet{Float64}}(undef, length(bixels))

        # Iterate through each control point in the field
        for beam in beams
            ΔMU = getΔMU(beam)
            gantry = getgantry(beam)
            mlc = getmlc(beam)
        
            beamlets .= Beamlet.(bixels, Ref(gantry))
        
            dose_fluence_matrix!(D, pos_fixed, beamlets, surf, calc)
            fluence!(Ψ, bixels, mlc)
        
            dose .+= ΔMU*D*vec(Ψ)
            
            if(show_progress)
                next!(p; showvalues=[(:MU, totalMU+=ΔMU)])
            end

        end
    end
    dose
end
