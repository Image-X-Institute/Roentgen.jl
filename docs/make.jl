using Documenter, DoseCalculations

makedocs(
    modules = [DoseCalculations],
    sitename="DoseCalculations.jl Documentation",
    pages = [
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",
        ]
    ]
    )

deploydocs(
    repo = "github.com/ACRF-Image-X-Institute/DoseCalculations.jl.git",
    )
