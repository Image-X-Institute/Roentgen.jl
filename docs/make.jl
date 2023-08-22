using Documenter, Roentgen

makedocs(
    modules = [Roentgen],
    sitename="Roentgen.jl Documentation",
    pages = [
        "index.md",
        "Dose Volume"=>[
            "DoseVolume.md",
            "DosePositions.md",
            "ExternalSurfaces.md"
        ],
        "Computing Dose"=>[
            "Fluence.md",
            "DoseFluenceMatrix.md"
        ],
        "Dose Calculation Algorithms"=>["ScaledIsoplaneKernel.md",],
        "Utilties"=>[
            "Structures.md"
        ],
        "API"=>["API.md",]
    ]
    )

deploydocs(
    repo = "github.com/Image-X-Institute/Roentgen.jl.git",
    )
