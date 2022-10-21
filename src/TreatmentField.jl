import Base.getindex, Base.lastindex, Base.length, Base.iterate, Base.show

export getmlc, getjaws, getϕg, getθb, getmeterset, getdoserate, getisocenter, getSAD, getmetersettype
export resample, getgantry, getΔMU
export VMATField

fixangle(angle) = angle > π ? angle - 2π : angle

"""
    Treatment Field

Abstract treatment field, basis for containing multiple treatment types (e.g. VMAT, IMRT, etc.)
"""
abstract type AbstractTreatmentField end

show_angle(angle, digits=2) = "$(round(rad2deg(angle); digits=digits))°"

"""
    iterate

Iteration of a treatment field, returning ControlPoint every time iterate is called.
"""
Base.iterate(field::AbstractTreatmentField) = Base.iterate(field, 1)
function Base.iterate(field::AbstractTreatmentField, i)
    if(length(field)<i)
        return nothing
    end
    field[i], i+1
end

Base.length(field::AbstractTreatmentField) = field.ncontrol

Base.eachindex(field::AbstractTreatmentField) = Base.OneTo(length(field))
Base.lastindex(field::AbstractTreatmentField) = length(field)

#--- Accessor methods ---------------------------------------------------------

getmlc(field::AbstractTreatmentField, args...) = field.mlc
getjaws(field::AbstractTreatmentField, args...) = field.jaws
getθb(field::AbstractTreatmentField, args...) = field.collimator_angle
getϕg(field::AbstractTreatmentField, args...) = field.gantry_angle
getSAD(field::AbstractTreatmentField, args...) = field.source_axis_distance
getmeterset(field::AbstractTreatmentField, args...) = field.meterset
getdoserate(field::AbstractTreatmentField, args...) = field.dose_rate
getisocenter(field::AbstractTreatmentField, args...) = field.isocenter

function Base.show(io::IO, field::AbstractTreatmentField)
    msg = "$(length(field)),"*
          " ϕg: $(show_angle(field.gantry_angle[1]))->$(show_angle(field.gantry_angle[end])),"*
          " θb: $(show_angle(field.collimator_angle)),"*
          " MU: $(field.meterset[end])"

    println(io, msg)
end

#--- Control Point ------------------------------------------------------------

struct ControlPoint{T<:AbstractFloat, TMLC<:AbstractMultiLeafCollimator} <: AbstractTreatmentField
    # Beam Limiting Devices
    mlc::TMLC
    jaws::Jaws{T} # Jaw Positions

    # Treatment Positions
    gantry_angle::T
    collimator_angle::T
    source_axis_distance::T

    ΔMU::T # Meterset   (MU)
    dose_rate::T        # Dose Rate (MU/s)
    isocenter::SVector{3, T}    # Isocenter Position
end

function Base.show(io::IO, pt::ControlPoint)
    msg = "ϕg: $(show_angle(getϕg(pt))),"*
          "θb: $(show_angle(getθb(pt))),"*
          "ΔMU: $(getΔMU(pt))"

    println(io, msg)
end

fixed_to_bld(pt::ControlPoint) = fixed_to_bld(getϕg(pt), getθb(pt), getSAD(pt))
getgantry(pt::ControlPoint) = GantryPosition(getϕg(pt), getθb(pt), getSAD(pt))
getΔMU(pt::ControlPoint) = pt.ΔMU

#--- VMAT --------------------------------------------------------------------------------------------------------------

abstract type AbstractVMATField <: AbstractTreatmentField end

struct VMATField{T<:AbstractFloat, TMLC} <: AbstractVMATField
    ncontrol::Int   # Number of Control Points

    # Beam Limiting Devices
    mlc::TMLC   # MLC Type

    jaws::Jaws{T}  # Jaw Positions

    # Treatment Positions
    gantry_angle::Vector{T}    # Gantry Angle (ʳ)
    collimator_angle::T  # Beam Limiting Device Angle (ʳ)
    source_axis_distance::T

    meterset::Vector{T}  # Cumulative Meterset   (MU)
    dose_rate::T   # Dose Rate (MU/s)

    isocenter::SVector{3, T} # Isocenter Position (in patient based coords)

    function VMATField(mlc, jaws, gantry_angles, collimator_angle, source_axis_distance, meterset, doserate, isocenter)
        ncontrol = length(gantry_angles)

        @assert length(mlc) == ncontrol
        @assert length(meterset) == ncontrol "length(meterset) ($(length(meterset))) != $ncontrol"
        new{typeof(doserate), typeof(mlc)}(ncontrol,
                                            mlc, jaws,
                                            gantry_angles, collimator_angle, source_axis_distance,
                                            meterset, doserate, isocenter)
    end
end

function Base.getindex(field::AbstractVMATField, i::Int)
    ControlPoint(getmlc(field, i),
                 getjaws(field),
                 getϕg(field, i),
                 getθb(field, i),
                 getSAD(field, i),
                 getΔMU(field, i),
                 getdoserate(field, i),
                 getisocenter(field, i))
end


function Base.getindex(field::AbstractVMATField, i::UnitRange{Int})
    VMATField(getmlc(field, i),
              getjaws(field),
              getϕg(field, i),
              getθb(field, i),
              getSAD(field, i),
              getmeterset(field, i),
              getdoserate(field, i),
              getisocenter(field, i))
end

Base.length(field::AbstractVMATField) = field.ncontrol

getgantry(field::AbstractVMATField) = GantryPosition.(getϕg(field),
                                                      getθb(field),
                                                      getSAD(field))

getgantry(field::AbstractVMATField, i::Int) = GantryPosition(getϕg(field, i),
                                                             getθb(field),
                                                             getSAD(field))

getmlc(field::AbstractVMATField, i) = field.mlc[i]
getϕg(field::AbstractVMATField, i) = field.gantry_angle[i]

getmeterset(field::AbstractVMATField, i) = field.meterset[i]

function getΔMU(field::AbstractVMATField, i::Int)
    i == 1 && return 0.5*(field.meterset[2]-field.meterset[1])
    i == length(field) && return 0.5*(field.meterset[end]-field.meterset[end-1])
    0.5*(field.meterset[i+1]-field.meterset[i-1])
end

#--- IO

"""
    save(file::HDF5.H5DataStore, field::VMATField)

Store `VMATField` data to an HDF5 file/group
"""
function save(file::HDF5.H5DataStore, field::VMATField)

    save(file, getmlc(field))
    save(file, getjaws(field))

    file["gantry_angle"] = getϕg(field)
    file["collimator_angle"] = getθb(field)
    file["source_axis_distance"] = getSAD(field)

    file["meterset"] = getmeterset(field)
    file["dose_rate"] = getdoserate(field)

    file["isocenter"] = Vector(getisocenter(field))

    nothing
end

"""
    load(file::HDF5.H5DataStore, field::VMATField)

Load `VMATField` data to an HDF5 file/group
"""
function load(::Type{VMATField}, file::HDF5.H5DataStore)
    mlc = load(MultiLeafCollimatorSequence, file)
    jaws = load(Jaws, file)

    names = ["gantry_angle", "collimator_angle", "source_axis_distance", "meterset", "dose_rate", "isocenter"]

    data = read.(getindex.(Ref(file), names))

    VMATField(mlc, jaws, data...)
end

#--- Resampling --------------------------------------------------------------------------------------------------------

"""
    resample(field, ηₛ::AbstractVector{T}; by=:time)

Resample a treatment field onto new times or meterset values.

Can resample either by time (`by=:time`) or MU (`by=:MU`, default).
"""
function resample(field::AbstractTreatmentField, ηₛ::AbstractVector{T}; by=:MU) where T<:AbstractFloat

    if(by==:time)
        η = field.meterset./field.dose_rate
    elseif(by==:MU)
        η = field.meterset
    end
    Δη = @. η[2:end]-η[1:end-1]

    ncontrol = length(ηₛ)
    ϕg = zeros(ncontrol)
    meterset = zeros(ncontrol)

    mlc = MultiLeafCollimatorSequence(ncontrol, getedges(field.mlc))

    for i in eachindex(ηₛ)
        j = min(length(Δη), searchsortedlast(η, ηₛ[i]))
        α = (ηₛ[i] - η[j])/Δη[j]

        # Interpolate Gantry Angle
        ϕg[i] = interp(field.gantry_angle[j], field.gantry_angle[j+1], α)
        meterset[i] = interp(field.meterset[j], field.meterset[j+1], α)

        # Interpolate MLC x
        mlc.positions[:, :, i] .= interp(field.mlc.positions[:,:,j], field.mlc.positions[:,:,j+1], α)        

    end
    VMATField(mlc, field.jaws,
              ϕg, field.collimator_angle, field.source_axis_distance,
              meterset, field.dose_rate,
              field.isocenter)
end

"""
    resample(field, Δη::T; by=:MU)

Resample at uniform steps `Δη`, from start to finish.

See `resample(field, ηₛ::AbstractVector{T}; by=:MU)` for details.
"""
function resample(field, Δη::T; by=:MU, include_end=true) where T<:Number

    if(by==:time)
        η_end = field.meterset[end]/field.dose_rate
    elseif(by==:MU)
        η_end = field.meterset[end]
    end

    η = collect(zero(T):Δη:η_end)
    if(include_end && η[end] < η_end)
        push!(η, η_end)
    end
    resample(field, η; by=by)
end
