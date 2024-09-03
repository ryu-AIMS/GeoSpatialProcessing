
"""
    to_multipolygon(raster::Raster{T, 2}) where {T<:Union{Bool,Int16}}

Convert raster to multipolygons.
Invert vertical axis.
"""
function to_multipolygon(
    raster::Raster{T, 2}
)::GI.Wrappers.MultiPolygon where {T<:Union{Bool,Int16}}
    return GO.polygonize(.==(0), raster[:, end:-1:1])
end

"""
    to_dataframe(mp::GI.Wrappers.MultiPolygon)::DataFrame

Create a DataFrame from multipolygons

"""
function to_dataframe(mp::GI.Wrappers.MultiPolygon)::DataFrame
    return DataFrame(geometry=mp.geom)
end

"""
    meters_to_degrees(x, lat)

Convert meters to degrees at a specified latitude value.
"""
function meters_to_degrees(x, lat)
    return x / (111.1*1000 * cosd(lat))
end

"""
    from_zero(v)

Ensures the coordinates (x and y) contained in vector `v` are scaled to start from 0.
"""
function from_zero(v)
    max_coord = maximum(v)
    if first(v) == max_coord
        v = reverse(v)
    end

    new_coords = [Vector{Union{Missing, Float64}}(missing, size(max_coord, 1)), Vector{Union{Missing, Float64}}(missing, size(max_coord, 1))]
    for (j, coords) in enumerate(new_coords)
        for val in eachindex(coords)
            coords[val] = max_coord[val] - v[j][val]
        end
    end

    return new_coords
end
