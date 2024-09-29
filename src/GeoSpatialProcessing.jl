module GeoSpatialProcessing

using Statistics

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
using GeometryBasics
using CoordinateTransformations

using Rasters
using DataFrames
import GeoDataFrames as GDF

using Distances: haversine
using
    Images,
    ImageFiltering,
    ImageMorphology
using
    GLMakie
    GeoMakie
    ColorSchemes

include("analysis/type_conversion.jl")
include("analysis/raster_processing.jl")
include("analysis/plotting.jl")

export
    to_multipolygon,
    to_dataframe
    set_consistent_missingval
    remove_orphaned_elements
    plot_map
    plot_map!
    plot_lines
end
