module SimplexGridFactory

using Printf: @printf, @sprintf
using LinearAlgebra: norm

using ElasticArrays: ElasticArray
import ExtendableGrids
using ExtendableGrids: dim_space, Coordinates, BFaceNodes, BFaceRegions, CellNodes, CellRegions
using GridVisualize: PyPlotType, plottertype, MakieType, GridVisualizer
using DocStringExtensions: SIGNATURES, TYPEDEF, TYPEDSIGNATURES
using FileIO: load
import MeshIO

include("utils.jl")

include("options.jl")

include("binnedpointlist.jl")
export BinnedPointList

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
export circle!, rect2d!, rect3d!, sphere!, bregions!, moveto!, lineto!
export model3d!, mesh3d!

include("plot.jl")
export builderplot

end # module
