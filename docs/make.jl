using Documenter, DoseCalculations

makedocs(
    modules = [DoseCalculations],
    sitename="DoseCalculations.jl Documentation",
    pages = [
        "Dose Volume"=>["ExternalSurfaces.md",
                        "Structures.md"],
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",],
        "API"=>["API.md",]
    ]
    )

deploydocs(
    repo = "github.com/ACRF-Image-X-Institute/DoseCalculations.jl.git",
    )
