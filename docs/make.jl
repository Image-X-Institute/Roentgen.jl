using Documenter
using DoseCalculations

makedocs(
    modules = [DoseCalculations],
    sitename="DoseCalculations.jl Documentation",
    pages = [
        "Dose Volume"=>["ExternalSurfaces.md",],
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",],
        "API"=>["API.md",]
    ]
    )
