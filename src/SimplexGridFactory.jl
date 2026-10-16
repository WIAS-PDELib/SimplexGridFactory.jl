"""
    SimplexGridFactory

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module SimplexGridFactory


using DocStringExtensions: SIGNATURES, TYPEDEF, TYPEDSIGNATURES
using ElasticArrays: ElasticArray
using ExtendableGrids: ExtendableGrids, dim_space, Coordinates, BFaceNodes, BFaceRegions, simplexgrid, BinnedPointList
using FileIO: load
using GridVisualize: GridVisualize, GridVisualizer, plot_triangulateio!
using LinearAlgebra: norm, det
using MeshIO: MeshIO
using Printf: @sprintf
using StaticArrays: @MArray

include("utils.jl")

include("options.jl")

include("simplexgridbuilder.jl")
export SimplexGridBuilder
export point!, facet!, polyfacet!, cellregion!, holepoint!, facetregion!, maxvolume!, regionpoint!, options!
export istriangulate, istetgen, maybewatertight
export flags

include("simplexgrid.jl")
export simplexgrid

include("triangle.jl")

include("tetgen.jl")

include("primitives.jl")
export circle!, arc!, rect2d!, rect3d!, sphere!, bregions!, moveto!, lineto!
export model3d!, mesh3d!

include("plot.jl")
export builderplot

end # module
