module DoseCalculations

using ProgressMeter
using Interpolations
import JSON

include("utils/FileConverters.jl")
include("utils/interpolation.jl")
include("utils/misc.jl")


include("CoordinateSystems.jl")

include("Gantry.jl")

include("BeamLimitingDevices/BeamLimitingDevices.jl")
include("DicomPlan.jl")

include("TreatmentField.jl")

include("Grids.jl")

include("ExternalSurfaces.jl")

include("DosePoints.jl")

include("Fluence.jl")
include("DoseCalculationAlgorithms/DoseCalculationAlgorithms.jl")

include("DoseFluenceMatrix.jl")

include("Structures.jl")
include("BeamLimitingDevicePlots.jl")

include("DoseReconstructions.jl")

end # module
