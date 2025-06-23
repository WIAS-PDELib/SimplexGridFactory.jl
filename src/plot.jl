"""
    builderplot(
        builder::SimplexGridBuilder;
        Plotter::Module = nothing,
        size = (650, 300),
        input_slot = (1, 1),
        output_slot = (1, 2),
        layout = (1, 2),
        vis = GridVisualizer(; Plotter, layout, size),
        circumcircles = false,
        reveal = true,
        kwargs...
    )

Two panel visualization of gridfactory with input and resulting grid
-  `builder` : Simplex grid builder
-  `Plotter`  : Plotter
-  `size`: size of plot
-  `input_slot`: slot in visualizer for input plot
-  `output_slot`: slot in visualizer for output plot
-  `layout`: layout of grid visualizer
-  `vis`: grid visualizer
-  `circumcircles`: plot circumcicles in output
-  `reveal`: reveal plot upon return
-  `kwargs...`: passed to output constructor; see [`default_options`](@ref) for available `kwargs`.
"""
builderplot(gb::SimplexGridBuilder; Plotter = nothing, kwargs...) = builderplot(gb, Plotter; kwargs...)

builderplot(builder::SimplexGridBuilder, ::Nothing; kwargs...) = nothing

function builderplot(
        builder::SimplexGridBuilder, Plotter::Module;
        size = (650, 300),
        input_slot = (1, 1),
        output_slot = (1, 2),
        layout = (1, 2),
        vis = GridVisualizer(; Plotter, layout, size),
        circumcircles = false,
        reveal = true,
        kwargs...
    )
    opts = blendoptions!(copy(builder.options); kwargs...)

    Triangulate = builder.Generator
    @assert(istriangulate(Triangulate))

    flags = makeflags(opts, :triangle)

    if opts[:verbose]
        @show flags
    end

    triin = nothing
    try
        triin = triangulateio(builder)
    catch err
        @error "Incomplete geometry description"
        rethrow(err)
    end

    if !isnothing(opts[:unsuitable])
        Triangulate.triunsuitable!(opts[:unsuitable])
    end

    triout, vorout = Triangulate.triangulate(flags, triin)

    plot_triangulateio!(
        vis[input_slot...],
        triin;
        title = "Input"
    )

    plot_triangulateio!(
        vis[output_slot...],
        triout;
        voronoi = length(vorout.pointlist) > 0 ? vorout : nothing,
        circumcircles,
        title = "Output"
    )

    if reveal
        return GridVisualize.reveal(vis)
    else
        return vis
    end
end
