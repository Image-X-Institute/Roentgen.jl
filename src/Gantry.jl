
struct GantryPosition{T}
    gantry_angle::T
    collimator_angle::T
    source_axis_distance::T
end

function Base.show(io::IO, gantry::GantryPosition)
    print(io, "ϕg=", show_angle(getϕg(gantry)), ", θb=", show_angle(getθb(gantry)), ", SAD=", getSAD(gantry))
end

getϕg(gantry::GantryPosition) = gantry.gantry_angle
getθb(gantry::GantryPosition) = gantry.collimator_angle
getSAD(gantry::GantryPosition) = gantry.source_axis_distance

fixed_to_bld(gantry::GantryPosition) = fixed_to_bld(getϕg(gantry), getθb(gantry), getSAD(gantry))
