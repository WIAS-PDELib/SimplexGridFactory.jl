
"""
$(TYPEDSIGNATURES)

Two panel visualization of gridfactory with input and resulting grid
See [`default_options`](@ref) for available `kwargs`.
"""
builderplot(gb::SimplexGridBuilder; Plotter = nothing, kwargs...) = builderplot(plottertype(Plotter), gb, Plotter; kwargs...)

builderplot(::Type{Nothing}, builder::SimplexGridBuilder, PyPlot; kwargs...) = nothing

function builderplot(::Type{PyPlotType}, builder::SimplexGridBuilder, PyPlot; kwargs...)
    p = GridVisualizer(; Plotter = PyPlot, layout = (1, 2), kwargs...)

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
        Triangulate.triunsuitable(opts[:unsuitable])
    end

    triout, vorout = Triangulate.triangulate(flags, triin)
    PyPlot.subplot(121)
    PyPlot.title("In")
    Triangulate.plot_triangulateio(PyPlot, triin)
    PyPlot.subplot(122)
    PyPlot.title("Out")
    Triangulate.plot_triangulateio(PyPlot, triout)
    PyPlot.tight_layout()
    PyPlot.gcf()
end
