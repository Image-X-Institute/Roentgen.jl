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
    snapped_range(x1, x2, Δ)

Create a range from x1 to x2 which is "snapped" to the step Δ.

Positions are "snapped" to the step value (e.g. a starting position of
x[1]-0.2Δx snaps to x[1]-Δx). The new range always includes the start and end
points of the original range

Examples:
- 0.1:1.:9.4 -> 0:1.:10.

"""
snapped_range(x1, x2, Δ) = Δ*(floor(Int, x1/Δ):ceil(Int, x2/Δ))
snapped_range(x::AbstractRange, Δ) = snapped_range(first(x), last(x), Δ)
