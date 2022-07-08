import Base.getindex, Base.lastindex, Base.length, Base.iterate, Base.show

export getmlcx, getmlc, getjaws, getϕg, getθb, getmeterset, getdoserate, getisocenter, getSAD, getmetersettype
export limit_to_jaws, resample, discretise

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

getmlcx(field::AbstractTreatmentField, args...) = field.mlcx
getmlc(field::AbstractTreatmentField, args...) = field.mlc
getjaws(field::AbstractTreatmentField, args...) = field.jaws
getθb(field::AbstractTreatmentField, args...) = field.θb
getϕg(field::AbstractTreatmentField, args...) = field.ϕg
getSAD(field::AbstractTreatmentField, args...) = field.SAD
getmeterset(field::AbstractTreatmentField, args...) = field.meterset
getmetersettype(field::AbstractTreatmentField, args...) = field.metersettype
getdoserate(field::AbstractTreatmentField, args...) = field.Ḋ
getisocenter(field::AbstractTreatmentField, args...) = field.isocenter

TreatmentField(mlcx, mlc, jaws, θb, ϕg::AbstractVector, meterset, Ḋ) = VMATField(mlcx, mlc, jaws, θb, ϕg, meterset, Ḋ)

struct ControlPoint{T<:AbstractFloat, TMLC<:MultiLeafCollimator, TArray<:AbstractArray{T}} <: AbstractTreatmentField
    mlcx::TArray
    mlc::TMLC

    jaws::Jaws{T} # Jaw Positions

    # Treatment Positions
    θb::T       # Beam Limiting Device Angle (ʳ)
    ϕg::T       # Gantry Angle (ʳ)
    meterset::T # Meterset   (MU)
    metersettype::Symbol   # Cumulative or discrete
    Ḋ::T        # Dose Rate (MU/s)
    isocenter::SVector{3, T}    # Isocenter Position
    SAD::T      # Source-Axis Distance (mm)

end

function Base.getindex(field::AbstractTreatmentField, index::Int) 
    ControlPoint(getmlcx(field, index),
                 getmlc(field),
                 getjaws(field),
                 getθb(field, index),
                 getϕg(field, index),
                 getmeterset(field, index),
                 getmetersettype(field),
                 getdoserate(field),
                 getisocenter(field),
                 getSAD(field))
end

#--- VMAT --------------------------------------------------------------------------------------------------------------

abstract type AbstractVMATField <: AbstractTreatmentField end

struct VMATField{T<:AbstractFloat, TMLC<:MultiLeafCollimator, TArray<:AbstractArray{T, 3}, TVec<:AbstractVector{T}} <: AbstractVMATField
    ncontrol::Int   # Number of Control Points

    # Beam Limiting Devices
    mlcx::TArray    # MLC Sequence
    mlc::TMLC   # MLC Type

    jaws::Jaws{T}  # Jaw Positions

    # Treatment Positions
    θb::T  # Beam Limiting Device Angle (ʳ)
    ϕg::TVec    # Gantry Angle (ʳ)
    meterset::TVec  # Cumulative Meterset   (MU)
    metersettype::Symbol   # Cumulative or discrete
    Ḋ::T   # Dose Rate (MU/s)

    isocenter::SVector{3, T} # Isocenter Position (in patient based coords)
    SAD::T
end

getmlcx(field::AbstractVMATField, i) = view(field.mlcx, :, :, i)

getθb(field::AbstractVMATField, i) = field.θb

getϕg(field::AbstractVMATField, i::AbstractVector) = view(field.ϕg, i)
getϕg(field::AbstractVMATField, i) = field.ϕg[i]

getmeterset(field::AbstractVMATField, i) = field.meterset[i]
getmeterset(field::AbstractVMATField, i::AbstractVector) = view(field.meterset, i)

function Base.getindex(field::AbstractVMATField, index::AbstractVector)
    VMATField(length(index),
              getmlcx(field, index),
              getmlc(field),
              getjaws(field),
              getθb(field),
              getϕg(field, index),
              getmeterset(field, index),
              getmetersettype(field),
              getdoserate(field),
              getisocenter(field),
              getSAD(field))
end

function Base.show(io::IO, field::AbstractTreatmentField)
    msg = "$(length(field)),"*
          " ϕg: $(show_angle(field.ϕg[1]))->$(show_angle(field.ϕg[end])),"*
          " θb: $(show_angle(field.θb)),"

    if(getmetersettype(field)==:cumulative)
        msg *= " MU: $(field.meterset[1])->$(field.meterset[end])"
    elseif(getmetersettype(field)==:discrete)
        msg *= " ΔMU: $(sum(field.meterset))"
    end

    println(io, msg)
end

Base.length(field::AbstractVMATField) = field.ncontrol


# Limit Plan to Jaws

"""
    limit_to_jaws(field::AbstractTreatmentField)

Apply to a `TreatmentField`, using the jaws as the lower and upper positions.
"""
function limit_to_jaws(field::AbstractTreatmentField)
    jaws = getjaws(field)
    mlc = getmlc(field)
    mlcx = getmlcx(field)

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
function resample(field, Δη::T; by=:time) where T<:Number

    if(by==:time)
        η_end = field.meterset[end]/field.Ḋ
    elseif(by==:MU)
        η_end = field.meterset[end]
    end

    η = collect(zero(T):Δη:η_end)
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
   
    mlcx = @. 0.5*(field.mlcx[:, :, 1:end-1] + field.mlcx[:, :, 2:end])
    ϕg = @. 0.5*(field.ϕg[1:end-1] + field.ϕg[2:end])

    ΔMU = field.meterset[2:end] - field.meterset[1:end-1]

    VMATField(ncontrol,
              mlcx, field.mlc, field.jaws,
              field.θb, ϕg,
              ΔMU, :discrete, field.Ḋ,
              field.isocenter, field.SAD)
end
