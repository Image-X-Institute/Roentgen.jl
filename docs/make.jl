using Documenter
using DoseCalculations

makedocs(
    modules = [DoseCalculations],
    sitename="DoseCalculations.jl Documentation",
    pages = [
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",
        ]
    ]
    )
