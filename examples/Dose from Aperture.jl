using Pkg; Pkg.activate("."); Pkg.instantiate()
using Roentgen
using SparseArrays

plan = load_dicom("path/to/dicom/RP.....dcm")
field = plan[1] # Select the first field
controlpoint = field[1] # Select the first control point

mesh = load_structure_from_ply("path/to/body.stl")
trans = patient_to_fixed(getisocenter(controlpoint))

surf = CylindricalSurface(transform(mesh, trans), 5.);

pos = DoseGridMasked(5., SurfaceBounds(surf), trans)
pos_fixed = trans.(pos);

calc = FinitePencilBeamKernel("path/to/kernel/file.jld")
calibrate!(calc, 100., 100., 1000.);

bixels = BixelGrid(getjaws(controlpoint), 5., 5.)
beamlets = Beamlet.(bixels, Ref(getgantry(controlpoint)));

Ψ = fluence(bixels, getmlc(controlpoint)); # Compute the fluence
ΔMU = getΔMU(controlpoint);

D = dose_fluence_matrix(SparseMatrixCSC, vec(pos), vec(beamlets), surf, calc);

dose = ΔMU*D*vec(Ψ);

# To VTK (Paraview)
write_vtk("dose", pos, "dose"=>dose)
# To NRRD (Slicer)
write_nrrd("dose.nrrd", pos, dose);
