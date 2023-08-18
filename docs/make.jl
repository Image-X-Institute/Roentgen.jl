using Documenter, Roentgen

makedocs(
    modules = [Roentgen],
    sitename="Roentgen.jl Documentation",
    pages = [
        "Dose Volume"=>[
            "DoseVolume.md",
            "DosePositions.md",
            "ExternalSurfaces.md",
            "Structures.md"],
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",],
        "API"=>["API.md",]
    ]
    )

deploydocs(
    repo = "github.com/ACRF-Image-X-Institute/Roentgen.jl.git",
    )
