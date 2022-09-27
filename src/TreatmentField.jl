import Base.getindex, Base.lastindex, Base.length, Base.iterate, Base.show

export getmlc, getjaws, getϕg, getθb, getmeterset, getdoserate, getisocenter, getSAD, getmetersettype
export limit_to_jaws, resample, discretise

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

Base.eachindex(field::AbstractTreatmentField) = Base.OneTo(length(field))
Base.lastindex(field::AbstractTreatmentField) = length(field)

#--- Accessor methods ---------------------------------------------------------

getmlc(field::AbstractTreatmentField, args...) = field.mlc
getjaws(field::AbstractTreatmentField, args...) = field.jaws
getθb(field::AbstractTreatmentField, args...) = field.collimator_angle
getϕg(field::AbstractTreatmentField, args...) = field.gantry_angle
getSAD(field::AbstractTreatmentField, args...) = field.source_axis_distance
getmeterset(field::AbstractTreatmentField, args...) = field.meterset
getmetersettype(field::AbstractTreatmentField, args...) = field.metersettype
getdoserate(field::AbstractTreatmentField, args...) = field.Ḋ
getisocenter(field::AbstractTreatmentField, args...) = field.isocenter

#--- Control Point ------------------------------------------------------------

struct ControlPoint{T<:AbstractFloat, TMLC<:AbstractMultiLeafCollimator} <: AbstractTreatmentField
    # Beam Limiting Devices
    mlc::TMLC
    jaws::Jaws{T} # Jaw Positions

    # Treatment Positions
    gantry_angle::T
    collimator_angle::T
    source_axis_distance::T

    meterset::T # Meterset   (MU)
    Ḋ::T        # Dose Rate (MU/s)
    isocenter::SVector{3, T}    # Isocenter Position
end

function Base.show(io::IO, pt::ControlPoint)
    msg = "ϕg: $(show_angle(getϕg(pt))),"*
          "θb: $(show_angle(getθb(pt))),"*
          "MU: $(getmeterset(pt))"

    println(io, msg)
end

fixed_to_bld(pt::ControlPoint) = fixed_to_bld(getϕg(pt), getθb(pt), getSAD(pt))

#--- VMAT --------------------------------------------------------------------------------------------------------------

abstract type AbstractVMATField <: AbstractTreatmentField end

struct VMATField{T<:AbstractFloat, TMLC<:AbstractMultiLeafCollimator} <: AbstractVMATField
    ncontrol::Int   # Number of Control Points

    # Beam Limiting Devices
    mlc::TMLC   # MLC Type

    jaws::Jaws{T}  # Jaw Positions

    # Treatment Positions
    gantry_angle::Vector{T}    # Gantry Angle (ʳ)
    collimator_angle::T  # Beam Limiting Device Angle (ʳ)
    source_axis_distance::T

    meterset::Vector{T}  # Cumulative Meterset   (MU)
    Ḋ::T   # Dose Rate (MU/s)

    isocenter::SVector{3, T} # Isocenter Position (in patient based coords)

    function VMATField(mlc, jaws, gantry_angles, collimator_angle, source_axis_distance, meterset, dose_rate, isocenter)
        ncontrol = length(gantry_angles)

        @assert length(mlc) == ncontrol
        @assert length(meterset) == ncontrol-1
        new{typeof(dose_rate), typeof(mlc)}(ncontrol,
                                            mlc, jaws,
                                            gantry_angles, collimator_angle, source_axis_distance,
                                            meterset, dose_rate, isocenter)
    end
end

function Base.getindex(field::VMATField, i::Int) 
    ControlPoint(field.mlc[i:i+1],
                 field.jaws,
                 0.5*(field.ϕg[i]+field.ϕg[i+1]),
                 field.θb,
                 field.SAD,
                 field.meterset[i],
                 field.Ḋ,
                 field.isocenter)
end

function Base.show(io::IO, field::AbstractTreatmentField)
    msg = "$(length(field)),"*
          " ϕg: $(show_angle(field.gantry_angle[1]))->$(show_angle(field.gantry_angle[end])),"*
          " θb: $(show_angle(field.collimator_angle)),"
          " MU: $(sum(field.meterset))"

    println(io, msg)
end

Base.length(field::AbstractVMATField) = field.ncontrol-1

getϕg(field::AbstractVMATField, i) = field.gantry_angle[i]
getmeterset(field::AbstractVMATField, i) = field.meterset[i]

# Limit Plan to Jaws

"""
    limit_to_jaws(field::AbstractTreatmentField)

Apply to a `TreatmentField`, using the jaws as the lower and upper positions.
"""
function limit_to_jaws(field::AbstractTreatmentField)
    jaws = getjaws(field)
    mlc = getmlc(field)

    i1, i2 = subset_indices(mlc, jaws.y)

    mlcx = mlcx[:, i1:i2, :]
    mlc = MultiLeafCollimator(mlc[i1:i2])

    VMATField(length(field), mlcx, mlc, jaws, field.θb, field.ϕg, field.meterset, field.metersettype, field.Ḋ, field.isocenter, field.SAD)
end

"""
    resample(field, ηₛ::AbstractVector{T}; by=:time)

Resample a treatment field onto new times or meterset values.

Can resample either by time (`by=:time`, default) or MU (`by=:MU`).
"""
function resample(field::AbstractTreatmentField, ηₛ::AbstractVector{T}; by=:time) where T<:AbstractFloat

    if(getmetersettype(field)!=:cumulative)
        throw(ErrorException("getmetersettype(field) != :cumulative"))
    end

    if(by==:time)
        η = field.meterset./field.Ḋ
    elseif(by==:MU)
        η = field.meterset
    end
    Δη = @. η[2:end]-η[1:end-1]

    ncontrol = length(ηₛ)
    mlcx = zeros(2, length(field.mlc), ncontrol)
    ϕg = zeros(ncontrol)
    meterset = zeros(ncontrol)

    mlcx = zeros(2, length(field.mlc), length(ηₛ))

    for i in eachindex(ηₛ)
        j = min(length(Δη), searchsortedlast(η, ηₛ[i]))
        α = (ηₛ[i] - η[j])/Δη[j]

        # Interpolate Gantry Angle
        ϕg[i] = interp(field.ϕg[j], field.ϕg[j+1], α)
        meterset[i] = interp(field.meterset[j], field.meterset[j+1], α)

        # Interpolate MLC x
        mlcxi = @view mlcx[:,:,i]
        mlcxi .= interp(field.mlcx[:,:,j], field.mlcx[:,:,j+1], α)        

        # for j=1:size(mlcxi, 2)
        #     if(mlcxi[1,j] <= field.jaws.x[1] || mlcxi[2,j] >= field.jaws.x[2])
        #         mlcxi[:, j] .= 0.
        #     end
        # end
    end

    # mlcx = @. max(field.jaws.x[1], min(field.jaws.x[2], mlcx))

    VMATField(ncontrol,
              mlcx, field.mlc, field.jaws,
              field.θb, ϕg,
              meterset, getmetersettype(field), field.Ḋ,
              field.isocenter, field.SAD)
end

"""
    resample(field, Δη::T; by=:time)

Resample at uniform steps `Δη`, from start to finish.

See `resample(field, ηₛ::AbstractVector{T}; by=:time)` for details.
"""
function resample(field, Δη::T; by=:time, include_end=true) where T<:Number

    if(by==:time)
        η_end = field.meterset[end]/field.Ḋ
    elseif(by==:MU)
        η_end = field.meterset[end]
    end

    η = collect(zero(T):Δη:η_end)
    if(include_end && η[end] < η_end)
        push!(η, η_end)
    end
    resample(field, η; by=by)
end


"""
    discretise(field)

Discretise a treatment field.

Returns a field of length `length(field)-1`. Two consecutive control points are
combined by averaging the apertures and differencing the meterset value.
"""
function discretise(field::AbstractTreatmentField)

    if(getmetersettype(field)!=:cumulative)
        throw(ErrorException("getmetersettype(field) != :cumulative"))
    end

    ncontrol = length(field)-1

    ϕg = @. 0.5*(field.ϕg[1:end-1] + field.ϕg[2:end])

    ΔMU = field.meterset[2:end] - field.meterset[1:end-1]

    VMATField(ncontrol,
              field.mlcx, field.mlc, field.jaws,
              field.θb, ϕg,
              ΔMU, :discrete, field.Ḋ,
              field.isocenter, field.SAD)
end
