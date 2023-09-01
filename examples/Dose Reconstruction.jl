using Pkg; Pkg.activate("."); Pkg.instantiate()
using Roentgen

plan = load_dicom("path/to/dicom/RP.....dcm")

mesh = load_structure_from_ply("path/to/body.stl")
trans = patient_to_fixed(getisocenter(plan[1]))

surf = CylindricalSurface(transform(mesh, trans), 5.);

pos = DoseGridMasked(5., SurfaceBounds(surf), trans);

vol = DoseVolume(pos, surf)

calc = FinitePencilBeamKernel("path/to/kernel/file.jld")
calibrate!(calc, 100., 100., 1000.);

dose = reconstruct_dose(vol, plan, calc; show_progress=false);

write_nrrd("dose.nrrd", pos, dose);
