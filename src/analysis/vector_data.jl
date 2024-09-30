using DataFrames, ArchGDAL, Statistics
import GeoDataFrames as GDF


"""
    find_intersections(
        x::DataFrame,
        y::DataFrame,
        x_id::Symbol,
        y_id::Symbol,
        units::String,
        threshold,
        y_geom_col::Symbol=:geometry;
        proportion::Bool=false,
        nearest::Bool=false
    )::DataFrame

Find the areas of `y` that intersect with each polygon in `x`.
`rel_areas` contains corresponding `y_id` for each intersecting polygon in x (can then be
joined to `x`).

If `proportion = true`: polygons of `y` are only chosen if the intersection with `x` is >
50% the area of `x`.

# Arguments
- `x` : The target GeoDataFrame to compare with
- `y` : GeoDataFrame containing polygons to match against `x`
- `xid` : Column name holding unique IDs for x geometries (referred to as GBRMPA_ID in rel_areas)
- `yid` : Column name holding variable of interest for y geometries
- `units` : String indicating the units to use for the nearest distance threshold
- `threshold` : Value to use for threshold for joining to nearest object when no overlap is found.
- `y_geom_col` : Column name holding geometries in y
- `proportion` : Only select y polygons if the intersection with x polygon is > 50% of x polygon area
                 (default: `false`).
- `nearest` : When there is no overlap join to nearest geometry instead of returning missing
                (within threshold), only specify units and threshold when nearest=true.
"""
function find_intersections(
    x::DataFrame,
    y::DataFrame,
    x_id::Symbol,
    y_id::Symbol,
    units::String,
    threshold,
    y_geom_col::Symbol=:geometry;
    proportion::Bool=false,
    nearest::Bool=false
)::DataFrame
    rel_areas = DataFrame(
        [Vector{Any}(missing, size(x, 1)) for _ in 1:2],
        [x_id, :area_ID]
    )
        for (x_i, reef_poly) in enumerate(eachrow(x))
            intersecting = DataFrame(
                [Vector{Any}(missing, size(y, 1)) for _ in 1:3],
                [x_id, :area_ID, :inter_area]
            )

            for (y_i, interest_area) in enumerate(eachrow(y))
                if ArchGDAL.intersects(reef_poly.geometry, interest_area[y_geom_col])
                    inter_area = ArchGDAL.intersection(
                        reef_poly.geometry, interest_area[y_geom_col]
                    )

                    inter_area = ArchGDAL.geomarea(inter_area)
                    if proportion
                        prop_area = inter_area / ArchGDAL.geomarea(reef_poly.geometry)

                        if prop_area >= 0.5
                            data = [reef_poly[x_id], interest_area[y_id], inter_area]

                        else
                            data = [missing, missing, missing]
                        end
                    else
                        data = [reef_poly[x_id], interest_area[y_id], inter_area]
                    end
                else
                    data = [reef_poly[x_id], missing, missing]
                end

                intersecting[y_i, :] = data
            end

            if all(ismissing, intersecting.area_ID)
                if nearest
                    distances = ArchGDAL.distance.([reef_poly.geometry], y[:, y_geom_col])

                    if units == "km"
                        distances .*= 111.1
                    end
                    if minimum(distances) < threshold
                        x_data = [intersecting[1, x_id], y[argmin(distances), y_id]]
                    else
                        x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
                    end
                else
                    x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
                end
            else
                dropmissing!(intersecting)
                max_inter_area = argmax(intersecting.inter_area)
                x_data = [intersecting[max_inter_area, x_id], intersecting[max_inter_area, :area_ID]]
            end

            rel_areas[x_i, :] = x_data
        end

    return rel_areas
end

"""
    rotate_geom(geom, degrees::Float64)

Rotate a polygon geometry by `degrees` in the clockwise direction.

# Arguments
- `geom` : Geometry to rotate. ArchGDAL or GeometryOPs origin.
- `degrees` : Angle in degrees to rotate by (clockwise rotation).
"""
function rotate_geom(geom, degrees::Float64)
    degrees == 0.0 && return geom

    theta = deg2rad(degrees)
    sinang, cosang = sincos(theta)

    # Center is used as pivot point
    cx, cy = GO.centroid(geom)

    # Extract points
    new_points = GI.coordinates(geom)
    #new_points = get_points(geom)

    rotate_point(p) = begin
        x, y = p
        x -= cx
        y -= cy
        new_x = x * cosang - y * sinang + cx
        new_y = x * sinang + y * cosang + cy
        SVector(new_x, new_y)
    end

    # Calculate new coordinates of each vertex
    @inbounds @simd for i in eachindex(new_points)
        new_points[i] = rotate_point(new_points[i])
    end

    return create_poly(new_points, GI.crs(geom))
end

"""
    move_geom(geom, new_centroid::Tuple)

Move a geom to a new centroid.

# Arguments
- `geom` : geometry to move
- `new_centroid` : Centroid given in lon, lat
"""
function move_geom(geom, new_centroid::Tuple)
    tf_lon, tf_lat = new_centroid .- GO.centroid(geom)
    f = CoordinateTransformations.Translation(tf_lon, tf_lat)
    return GO.transform(f, geom)
end

"""
    find_horiz(geom)

Locate the first horizontal line contained in `geom`.
"""
function find_horiz(geom)
    coords = collect(GI.coordinates(geom)...)
    first_coord = first(coords)
    second_coord = coords[
        (getindex.(coords, 2) .∈ first_coord[2]) .&&
        (getindex.(coords, 1) .∉ first_coord[1])
    ]

    return [tuple(first_coord...), tuple(first(second_coord)...)]
end

"""
    angle_cust(a, b)

Find the angle between vectors `a` and `b`.
When using for geospatial processing it is recommended to use `from_zero` function first,
particularly when using a meters projection system.
"""
function angle_cust(a, b)
    return acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

"""
    simplify_exclusions!(exclusions::DataFrame, simplify_tol::Real)::DataFrame

Simplify exclusions by applying convex hull, and simplification operations to polygons.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons. Default = 0.1
"""
function simplify_exclusions!(exclusions::DataFrame; simplify_tol::Real)
    exclusions.geometry .= AG.convexhull.(exclusions.geometry)
    exclusions.geometry .= AG.simplify.(exclusions.geometry, simplify_tol)
    return exclusions
end

"""
    buffer_exclusions!(exclusions::DataFrame, buffer_dist::Real)::DataFrame

Buffer exclusion zones by a specified distance.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `buffer_dist::Real`: The buffer distance. Default = 1.0
"""
function buffer_exclusions!(exclusions::DataFrame; buffer_dist::Real)
    exclusions.geometry .= AG.buffer.(exclusions.geometry, buffer_dist)
    return exclusions
end

"""
    unionize_overlaps!(exclusions::DataFrame)::DataFrame

Unionize overlapping exclusion zones.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
"""
function unionize_overlaps!(exclusions::DataFrame)
    geometries = exclusions.geometry
    n = length(geometries)

    for i in 1:n
        geom1 = geometries[i]

        for j in i+1:n
            geom2 = geometries[j]

            if AG.overlaps(geom1, geom2)
                union_geom = AG.union(geom1, geom2)
                exclusions.geometry[i] = union_geom
                exclusions.geometry[j] = union_geom

                for k in 1:n
                    if AG.overlaps(union_geom, geometries[k])
                        exclusions.geometry[k] = union_geom
                    end
                end
            end

            if AG.contains(geom1, geom2)
                exclusions.geometry[j] = geom1
            elseif AG.contains(geom2, geom1)
                exclusions.geometry[i] = geom2
            end
        end
    end

    # Remove duplicate unionized geometries
    unique_geometries = unique(exclusions.geometry[.!AG.isempty.(exclusions.geometry)])
    exclusions = DataFrame(geometry = unique_geometries)

    return exclusions
end

"""
    polygon_to_lines(
        polygon::Union{Vector{T},T,GIWrap.MultiPolygon}
    ) where {T<:GIWrap.Polygon}

Extract the individual lines between vertices that make up the outline of a polygon.
"""
function polygon_to_lines(
    polygon::Union{Vector{T},T,GIWrap.MultiPolygon}
) where {T<:GIWrap.Polygon}
    poly_lines = [
        GO.LineString(GO.Point.(vcat(GI.getpoint(geometry)...)))
        for geometry in polygon.geom
    ]

    return vcat(poly_lines...)
end
