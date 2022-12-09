export GantryPosition

"""
    GantryPosition{T}

The position/rotation of the gantry and beam-limiting device.

Stores:
- `gantry_angle`: As defined by IEC
- `collimator_angle`: As defined by IEC (beam limiting device angle)
- `source_axis_distance`: Distance between source and isocenter
- `position`: Position of the source in the IEC Fixed coordinates
"""
struct GantryPosition{T}
    gantry_angle::T
    collimator_angle::T
    source_axis_distance::T
    position::SVector{3, T}

    function GantryPosition(ϕg, θb, SAD)
        ϕg, θb, SAD = promote(ϕg, θb, SAD)
        T = typeof(ϕg)
        pos = SAD*SVector(sin(ϕg), zero(T), cos(ϕg))
        new{T}(ϕg, θb, SAD, pos)
    end
end

function Base.show(io::IO, gantry::GantryPosition)
    print(io, "ϕg=", show_angle(getϕg(gantry)), ", θb=", show_angle(getθb(gantry)), ", SAD=", getSAD(gantry))
end

getϕg(gantry::GantryPosition) = gantry.gantry_angle
getθb(gantry::GantryPosition) = gantry.collimator_angle
getSAD(gantry::GantryPosition) = gantry.source_axis_distance

getposition(gantry::GantryPosition) = gantry.position

fixed_to_bld(gantry::GantryPosition) = fixed_to_bld(getϕg(gantry), getθb(gantry), getSAD(gantry))
bld_to_fixed(gantry::GantryPosition) = bld_to_fixed(getϕg(gantry), getθb(gantry), getSAD(gantry))
