


export show_supertypes

show_supertypes(x) = show_supertypes(typeof(x))

function show_supertypes(T::DataType)
    type_list = DataType[T]
    while T != Any
        T = supertype(T)
        push!(type_list, T)
    end
    for (i, T) in enumerate(Iterators.reverse(type_list))
        println(" "^max(0, i-2) * "→ "^min(i-1, 1), T)
    end
end

"""
    scale_to_isoplane(pᵢ, z_plane)

Scale the position `pᵢ` to the plane at distance `z_plane`.
"""
function scale_to_isoplane(pᵢ, z_plane)
    x, y, z = pᵢ
    α = z_plane/z
    SVector(α*x, α*y)
end

"""
    minmax(x, l, u)

Return `x` if `l<=x<=u`, `l` if `x<l` or `u` if `u<x`
"""
minmax(x, l, u) = max(l, min(x, u)) # Ensures l<=x<=u