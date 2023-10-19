"""
    Dose Reconstructions

High-level functions to compute dose reconstructions from a given set of beams.
"""

export reconstruct_dose

"""
    reconstruct_dose(vol::AbstractDoseVolume, plan::AbstractTreatmentPlan,
        calc::AbstractDoseAlgorithm; Δx=5., show_progress=true, maxradius=100.)

Dose reconstruction from a treatment plan.

Requires a dose volume (`vol`), a treatment plan (`plan`), and a dose calculation
algorithm (`calc`).

Optional arguments include:
- `Δx`: Size of each bixel in the fluence grid (defaults to 5.)
- `show_progress`: If true (default), displays the progress
- `maxradius`: Kernel truncation radius
"""
function reconstruct_dose(vol::AbstractDoseVolume, plan::AbstractTreatmentPlan,
    calc::AbstractDoseAlgorithm; Δx=5., show_progress=true, maxradius=100.)

    pos = getpositions(vol)
    surf = getsurface(vol)

    dose = zeros(length(pos))
    Ψ = zeros(1)

    if(show_progress)
        totalMU = 0.
        progress = Progress(sum(length.(plan)))
    end

    # Iterate through each field in the plan
    for field in plan

        pos_fixed = patient_to_fixed(getisocenter(field)).(pos)

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
                                                 vec(beamlets), surf, calc;
                                                 maxradius=maxradius);
                             dims=2)
            
            if(show_progress)
                next!(progress; showvalues=[(:MU, totalMU+=ΔMU)])
            end
        end

    end
    dose
end
