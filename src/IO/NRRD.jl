
function _get_nrrd_props(pos::AbstractDoseGrid)
    x, y, z = getaxes(pos)
    origin = x[1], y[1], z[1]
    directions = (step(x), 0, 0), (0, step(y), 0), (0, 0, step(z))
    
    Dict("space origin"=>origin,
         "space directions"=>directions,
         "space"=>"left-posterior-superior",
         "kinds"=>("domain", "domain", "domain"))
end

"""
    write_nrrd(filename::String, pos::DoseGrid, data::AbstractArray)

Write `DoseGrid` to the NRRD data format.
"""
function write_nrrd(filename::String, pos::DoseGrid, data::AbstractArray)
    props = _get_nrrd_props(pos::AbstractDoseGrid)
    FileIO.save(filename, reshape(data, size(pos)), props=props)
end

"""
    write_nrrd(filename::String, pos::DoseGridMasked, data::AbstractArray; fillvalue=NaN)

Write `DoseGridMasked` to the NRRD data format.

Masked points in `DoseGridMasked` are filled with `fillvalue`, defaults to `NaN`.
"""
function write_nrrd(filename::String, pos::DoseGridMasked, data::AbstractArray; fillvalue=NaN)

    props = _get_nrrd_props(pos)

    datagrid = fill(fillvalue, gridsize(pos)...)
    for i in eachindex(data)
        I = pos.indices[i]
        datagrid[I] = data[i]
    end
    
    FileIO.save(filename, datagrid, props=props)
end
