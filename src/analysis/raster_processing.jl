
"""
    set_consistent_missingval(raster, val)

Replace value used to indicate no data, and return a Raster type with this value set.

## Note
`raster` should cover the smallest possible extent to reduce processing time.

# Arguments
- `raster` : data to reconfigure
- `val` : Value to use to indicate no data

# Returns
Raster
"""
function set_consistent_missingval(raster, val)
    replace_missing!(raster, val)
    return Raster(raster; missingval=val)
end

"""
    remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})

Cleans up valid pixels that are by themselves and not worth including in later assessments.

# Arguments
- `rst_mask` : Mask of valid raster locations
- `min_cluster_size` : Number of elements that need to be clustered together to be kept
- `box_size` : area to search around center pixel (width, height). Must be odd numbers.
"""
function remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})
    labels = label_components(rst_mask, strel_box(box_size))

    # Count the size of each component
    component_sizes = component_lengths(labels)

    # Mask components to keep
    keep_mask = component_sizes .>= min_cluster_size

    cleaned_raster = map(x -> keep_mask[x], labels) .* rst_mask

    return cleaned_raster
end

"""
    filter_distances(
        target_rast::Raster,
        dist_buffer
    )::Raster

Apply a mask to exclude pixels that are outside the indicated distance buffer(s).

`target_rast` and the `dist_buffer` should be in the same CRS (e.g., EPSG:7844 / GDA2020).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
- `dist_buffer` : Buffer geometry to use as the mask.

# Returns
- Masked boolean raster indicating pixels that are within the target distance.
"""
function filter_distances(target_rast::Raster, dist_buffer)::Raster
    # Mask out areas outside considered distance from port
    return mask(Raster(target_rast; missingval=0); with=dist_buffer)
end

"""
    calc_distances(
        target_rast::Raster,
        gdf::DataFrame;
        units::String="NM"
    )::Raster

Calculate the minimum distance from each point in gdf.geometry to each valid pixel in target_rast.

`target_rast` and the `gdf` should be in the same CRS (e.g., EPSG:7844 / GDA2020).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to calculate distances from.
- `gdf` : GeoDataFrame of 'points' for distance calculation.

# Returns
- Raster of distances from each cell to the closest point in gdf.
"""
function calc_distances(
    target_rast::Raster,
    gdf::DataFrame;
    units::String="NM"
)::Raster
    tmp_areas = Float32.(copy(target_rast))

    # First dimension is the rows (longitude)
    # Second dimension is the cols (latitude)
    raster_lon = Vector{Float64}(tmp_areas.dims[1].val)
    raster_lat = Vector{Float64}(tmp_areas.dims[2].val)

    @floop for (lon_ind, lon) in enumerate(raster_lon)
        for (lat_ind, lat) in enumerate(raster_lat)
            if tmp_areas[lon_ind, lat_ind] != 0.0
                point = AG.createpoint()
                AG.addpoint!(point, lon, lat)

                pixel_dists = AG.distance.([point], gdf.geometry)
                geom_point = gdf[argmin(pixel_dists), :geometry]
                geom_point = (AG.getx(geom_point, 0), AG.gety(geom_point, 0))

                dist_nearest = Distances.haversine(geom_point, (lon, lat))

                # Convert from meters to nautical miles
                if units == "NM"
                    dist_nearest = dist_nearest / 1852
                end

                # Convert from meters to kilometers
                if units == "km"
                    dist_nearest = dist_nearest / 1000
                end

                tmp_areas.data[lon_ind, lat_ind] = Float32(dist_nearest)
            end
        end
    end

    tmp_areas = rebuild(tmp_areas, missingval=Float32(0.0))
    return tmp_areas
end

"""
    resample_and_write(
        input_raster::Union{Raster,Nothing},
        template_raster::Raster
        dst_file::String
    )::Nothing

Resample `input_raster` to `template_raster` to ensure matching spatial extent, CRS and resolution.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `input_raster` : Input raster dataset for resampling to template_raster.
- `rst_template` : Template raster for resampling.
- `dst_file` : File location name to create output file. Should include variable and region information.
"""
function resample_and_write(
    input_raster::Union{Raster,Nothing},
    rst_template::Raster,
    dst_file::String
)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    resample(input_raster; to=rst_template, filename=dst_file, format="COG")

    return nothing
end
