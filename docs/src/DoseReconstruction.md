
# Dose Reconstruction from a Treatment Plan

Roentgen.jl provides a high-level interface to reconstructing dose, enabling quick and easy dose reconstructions.

[`reconstruct_dose`](@ref) takes a [dose volume](@ref DoseVolumes), [Treatment Plan](@ref) and a [dose calculation algorithm](@ref DoseAlgorithms), and computes the dose for all fields in the plan.
An example of this is provided in the [Dose Reconstruction](https://github.com/Image-X-Institute/Roentgen.jl/blob/main/examples/Dose%20Reconstruction.ipynb) notebook.

```@docs
reconstruct_dose
```