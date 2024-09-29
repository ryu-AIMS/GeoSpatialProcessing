function _convert_plottable(gdf::Union{DataFrame,DataFrameRow}, geom_col::Symbol)
    local plottable
    try
        if gdf isa DataFrame
            plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))
        else
            plottable = GeoMakie.geo2basic(AG.forceto(gdf[geom_col], AG.wkbPolygon))
        end
    catch
        # Column is already in a plottable form, or some unrelated error occurred
        if gdf isa DataFrame
            plottable = gdf[:, geom_col]
        else
            plottable = [gdf[geom_col]]
        end
    end

    return plottable
end

"""
    plot_map(gdf::DataFrame; geom_col::Symbol=:geometry, color_by::Symbol)

Convenience plot function.

# Arguments
- `gdf` : GeoDataFrame
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `geom_col` : Column name holding geometries to plot
"""
function plot_map(gdf::Union{DataFrame,DataFrameRow}; geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)
    poly!(ga, plottable)

    display(f)

    return f, ga
end

function plot_map!(ga::GeoAxis, gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing

    plottable = _convert_plottable(gdf, geom_col)
    if !isnothing(color)
        poly!(ga, plottable, color=color)
    else
        poly!(ga, plottable)
    end

    # Set figure limits explicitly
    xlims!(ga)
    ylims!(ga)

    return nothing
end

function plot_map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
    return plot_map!(current_axis(), gdf; geom_col=geom_col, color=color)
end

function colorscheme_alpha(cscheme::ColorScheme, alpha = 0.5)
    ncolors = length(cscheme)

    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

function plot_map(
    gdf::Union{DataFrame,DataFrameRow},
    color_by::Symbol,
    alpha=1;
    geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(df[:, color_by]),1) == 2
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)[2:3]
        palette = ColorSchemes.tableau_20.colors[2:3]
    elseif size(unique(df[:, color_by]),1) <= 20
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
        palette = ColorSchemes.tableau_20.colors
    else
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    polys = poly!(ga, plottable, color=palette[color_indices])

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
    return f, ga
end
