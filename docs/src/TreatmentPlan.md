# Treatment Plan

Roentgen.jl supports the usage of treatment plans.
These are typically stored in the DICOM RP file format.

Treatment plans consist of a number of treatment fields, which specify machine patient setup parameters.

Each treatment field is essentially a vector of control points (*i.e.* `Vector{ControlPoint}`).
Indexing a treatment field will return a [`Roentgen.ControlPoint`](@ref),
```julia
control_point = field[1]
```
Iterating a field will return each control point sequentially.

## DICOM

Treatment plans can be loaded in from a DICOM RP file using,
```julia
plan = load_dicom("path/to/RP....dcm")
```
