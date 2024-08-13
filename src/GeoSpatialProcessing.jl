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


include("analysis/type_conversion.jl")

export
    to_multipolygon,
    to_dataframe


end
