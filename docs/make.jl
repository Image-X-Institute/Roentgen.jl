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
            "DoseReconstruction.md",
            "Fluence.md",
            "DoseFluenceMatrix.md",
            "DoseCalculationAlgorithms.md",
            ],
        "Utilties"=>[
            "Structures.md",
            "TreatmentPlan.md"
        ],
        "API.md"
    ]
    )

deploydocs(
    repo = "github.com/Image-X-Institute/Roentgen.jl.git",
    )
