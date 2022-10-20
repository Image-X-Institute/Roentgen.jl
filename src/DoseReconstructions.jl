"""
    Dose Reconstructions

High-level functions to compute dose reconstructions from a given set of beams.
"""

export reconstruct_dose

"""
    dose_from_beam(pos, surf, beam, calc)

Compute the dose from a beam.
"""
function dose_from_beam(pos, surf, beam, calc, bixels, D, Ψ)
    ΔMU = getΔMU(beam)

    mlc = getmlc(beam)

    dose_fluence_matrix!(D, pos, bixels, getgantry(beam), surf, calc)
    fluence!(Ψ, bixels, mlc)

    ΔMU*D'*Ψ
end


"""
    reconstruct_dose(pos, surf, plan, calc; Δx=1., ΔMU=2., show_progess=true)

Reconstruct the dose.

Requires dose positions (`pos`), an external surface (`surf`), a set of beams
(`plan`), and the dose calculation algorithm (`calc`). Optional arguments
include:
- `Δx`: Size of each bixel in the fluence grid (defaults to 1.)
- `ΔMU`: Meterset increments (defaults to 2.)
- `show_progess`: Whether to display the progress (defaults to `true`)
"""
function reconstruct_dose(pos, surf, plan, calc; Δx=1., ΔMU=2., show_progess=true)
    # Allocate dose arrays
    dose = zeros(length(pos))#, Threads.nthreads())


    # Iterate through each field in the plan
    for field in plan

        beams = resample(field, ΔMU, by=:MU)
        
        bixels = vec(bixel_grid(getmlc(field, 1), getjaws(field), Δx))
        D = spzeros(length(bixels), length(pos))
        Ψ = zeros(size(bixels))

        pos_fixed = patient_to_fixed(getisocenter(field)).(pos)

        if(show_progess)
            totalMU = 0.
            p = Progress(length(beams))
        end

        # Iterate through each control point in the field
        for i in eachindex(beams)
            dose .+= dose_from_beam(pos_fixed, surf, beams[i], calc, bixels, D, Ψ)
            # break
            
            if(show_progess)
                next!(p; showvalues=[(:MU, totalMU+=getΔMU(beams, i))])
            end

        end
    end
    dose
end
