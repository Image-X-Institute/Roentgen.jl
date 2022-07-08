"""
    Dose Reconstructions

High-level functions to compute dose reconstructions from a given set of beams.
"""

export reconstruct_dose

"""
    dose_from_beam(pos, surf, beam, calc)

Compute the dose from a beam.
"""
function dose_from_beam(pos, surf, beam, calc)
    ΔMU = getmeterset(beam)

    bixels = bixels_from_bld(getmlcx(beam), getmlc(beam), getjaws(beam))

    trans = fixed_to_bld(getϕg(beam), getθb(beam), getSAD(beam))
    pos_bld = [trans(pos[i]) for i in eachindex(pos)]

    D = dose_fluence_matrix(pos_bld, bixels, surf, calc)

    ΔMU*vec(sum(D, dims=1))
end


"""
    reconstruct_dose(pos, surf, plan, calc; bixel_size=(1., 5.), ΔMU=2., show_progess=true)

Reconstruct the dose.

Requires dose positions (`pos`), an external surface (`surf`), a set of beams
(`plan`), and the dose calculation algorithm (`calc`). Optional arguments
include:
- `bixel_size`: Size of each bixel in the fluence grid (defaults to 1, 5)
- `ΔMU`: Meterset increments
- `show_progess`: Whether to display the progress (defaults to `true`)
"""
function reconstruct_dose(pos, surf, plan, calc; bixel_size=(1., 5.), ΔMU=2., show_progess=true)
    # Allocate dose arrays
    dose = zeros(length(pos))#, Threads.nthreads())

    # Iterate through each field in the plan
    for (n, field) in enumerate(plan)

        beams = discretise(resample(field, ΔMU, by=:MU))

        if(show_progess) p = Progress(length(beams)) end

        # Iterate through each control point in the field
        for i in eachindex(beams) 
            dose .+= dose_from_beam(pos, surf, beams[i], calc)
            
            if(show_progess) next!(p) end

        end
    end
    vec(sum(dose, dims=2))
end
