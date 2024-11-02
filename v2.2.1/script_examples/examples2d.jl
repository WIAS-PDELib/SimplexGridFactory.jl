#
# 2D Triangulations
# =================
#
# These examples can be loaded into Julia (Revise.jl recommended)
#
# These examples return either a [`SimplexGridBuilder`](@ref) struct which can
# be turned into a grid by calling [`simplexgrid`](@ref)(builder), or
# just an ExtendableGrid with default simplex grid data.
#
# The control flags for [Triangle](https://github.com/JuliaGeometry/Triangulate.jl)
# are created based on default options
# provided by this module which try to ensure "good" grids for
# FEM and FVM computations. These are documented in [`SimplexGridFactory.default_options`](@ref)
# Occasional [`options!`](@ref) statements in the examples overwrite these defaults.
#
#
# This test code is released under the license conditions of
# Triangulate.jl
#

using SimplexGridFactory
using ExtendableGrids
using LinearAlgebra
using Triangulate

# ## Domain triangulation
# Here we just describe a domain as a polygon and
# mesh it.
#
function triangulation_of_domain()
    builder = SimplexGridBuilder(; Generator = Triangulate)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 1, 0)
    p3 = point!(builder, 1, 2)
    p4 = point!(builder, 0, 1)
    p5 = point!(builder, -1, 2)

    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facet!(builder, p3, p4)
    facet!(builder, p4, p5)
    facet!(builder, p5, p1)

    builder
end

#
# ![](triangulation_of_domain.png)
#

#
# ## Triangulation with size control and boundary markers
# The previous example was a little bit bland.
# We miss:
# - size control for the triangles
# - differently marking of boundary parts
# 
function nicer_triangulation_of_domain()
    builder = SimplexGridBuilder(; Generator = Triangulate)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 1, 0)
    p3 = point!(builder, 1, 2)
    p4 = point!(builder, 0, 1)
    p5 = point!(builder, -1, 2)

    facetregion!(builder, 1)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facetregion!(builder, 2)
    facet!(builder, p3, p4)
    facet!(builder, p4, p5)
    facetregion!(builder, 3)
    facet!(builder, p5, p1)

    options!(builder; maxvolume = 0.01)

    builder
end
#
# ![](nicer_triangulation_of_domain.png)
#

# ## Triangulation with subregions
# Here we create different subregions and
# apply the maxvolume constraint to
# the subregions
#
function triangulation_of_domain_with_subregions()
    builder = SimplexGridBuilder(; Generator = Triangulate)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 1, 0)
    p3 = point!(builder, 1, 2)
    p4 = point!(builder, 0, 1)
    p5 = point!(builder, -1, 2)

    facetregion!(builder, 1)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facetregion!(builder, 2)
    facet!(builder, p3, p4)
    facet!(builder, p4, p5)
    facetregion!(builder, 3)
    facet!(builder, p5, p1)

    facetregion!(builder, 4)
    facet!(builder, p1, p4)

    cellregion!(builder, 2)
    maxvolume!(builder, 0.1)
    regionpoint!(builder, -0.1, 0.5)

    cellregion!(builder, 3)
    maxvolume!(builder, 0.01)
    regionpoint!(builder, 0.2, 0.2)

    builder
end
#
# ![](triangulation_of_domain_with_subregions.png)
#

# ## Direct specification of input arrays
# Of course we can specify the input for Triangle
# directly. The aim of SimplexBuilder is to avoid
# the tedious and error prone counting connected
# with this approach.
function direct_square(Generator = Triangulate)
    simplexgrid(Generator;
                points = [0 0; 0 1; 1 1; 1 0]',
                bfaces = [1 2; 2 3; 3 4; 4 1]',
                bfaceregions = [1, 2, 3, 4],
                regionpoints = [0.5 0.5;]',
                regionnumbers = [1],
                regionvolumes = [0.01])
end
#
# ![](direct_square.png)
#

# ## Local refinement
#
# We can interface to Triangle's unsuitable mechanism
#
function square_localref()
    builder = SimplexGridBuilder(; Generator = Triangulate)
    cellregion!(builder, 1)
    maxvolume!(builder, 0.01)
    regionpoint!(builder, 0.5, 0.5)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 1, 0)
    p3 = point!(builder, 1, 1)
    p4 = point!(builder, 0, 1)

    facetregion!(builder, 1)
    facet!(builder, p1, p2)
    facetregion!(builder, 2)
    facet!(builder, p2, p3)
    facetregion!(builder, 3)
    facet!(builder, p3, p4)
    facetregion!(builder, 4)
    facet!(builder, p4, p1)

    refinement_center = [0.5, 0.5]
    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y2 + y2 + y3) / 3]
        dist = norm(bary - refinement_center)
        if area > 0.01 * dist
            return 1
        else
            return 0
        end
    end
    options!(builder; unsuitable = unsuitable)
    builder
end
#
# ![](square_localref.png)
#

# ## Domain with holes
# We can generate domains with holes.
# This at once shall demonstrate how the chosen 
# API approach eases bookeeping of features added to the
# geometry description
#
function swiss_cheese_2d()
    function circlehole!(builder, center, radius; n = 20)
        points = [point!(builder, center[1] + radius * sin(t), center[2] + radius * cos(t)) for t in range(0, 2π; length = n)]
        for i = 1:(n - 1)
            facet!(builder, points[i], points[i + 1])
        end
        facet!(builder, points[end], points[1])
        holepoint!(builder, center)
    end

    builder = SimplexGridBuilder(; Generator = Triangulate)
    cellregion!(builder, 1)
    maxvolume!(builder, 0.1)
    regionpoint!(builder, 0.1, 0.1)

    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 10, 0)
    p3 = point!(builder, 10, 10)
    p4 = point!(builder, 0, 10)

    facetregion!(builder, 1)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facet!(builder, p3, p4)
    facet!(builder, p4, p1)

    holes = [8.0 4.0;
             1.0 2.0;
             8.0 9.0;
             3.0 4.0;
             4.0 6.0;
             7.0 9.0;
             4.0 7.0;
             7.0 5.0;
             2.0 1.0;
             4.0 1.0;
             4.0 8.0;
             2.0 8.0;
             3.0 6.0;
             4.0 9.0;
             9.0 1.0;
             9.0 1.0;
             6.0 9.0;
             8.0 9.0;
             3.0 5.0;
             1.0 4.0]'

    radii = [0.15, 0.15, 0.1, 0.35, 0.2, 0.3, 0.1, 0.4, 0.1, 0.4, 0.4, 0.15, 0.2, 0.2, 0.2, 0.35, 0.15, 0.25, 0.15, 0.25]

    for i = 1:length(radii)
        facetregion!(builder, i + 1)
        circlehole!(builder, holes[:, i], radii[i])
    end

    builder
end
#
# ![](swiss_cheese_2d.png)
#


# ## Remeshing another grid
#
# The `bregions!` method allows to use another grid as geometry description
# 
function remesh_2d()
    b = SimplexGridBuilder(; Generator = Triangulate)
    X=0:0.1:1
    grid1 = simplexgrid(X, X)
    bregions!(b,grid1)
    simplexgrid(b,maxvolume=0.01)
end
#
# ![](remesh_2d.png)
#

# ## Glueing in another grid
#
# The `bregions!` method allows to extract parts of the geometry description from
# an already existing grid.
# 
function glue_2d()
    b = SimplexGridBuilder(; Generator = Triangulate)

    h = 0.25
    X = 0:h:10
    Y = 0:h:4

    grid1 = simplexgrid(X, Y)
    bfacemask!(grid1, [0, 2], [10, 2], 7)
    cellmask!(grid1, [0, 0], [10, 2], 2)

    xmin = -2
    x0min = 0
    x0max = 10
    xmax = 12
    y0min = 0
    ymin = -2
    y0max = 4
    ymax = 6
    facetregion!(b, 8)
    facet!(b, (xmin, ymin), (xmax, ymin))
    facet!(b, (xmax, ymin), (xmax, 2))
    facet!(b, (xmax, 2), (xmax, ymax))
    facet!(b, (xmin, ymax), (xmax, ymax))
    facet!(b, (xmin, ymin), (xmin, 2))
    facet!(b, (xmin, 2), (xmin, ymax))

    facetregion!(b, 7)
    facet!(b, (x0min, 2), (xmin, 2))
    facet!(b, (x0max, 2), (xmax, 2))

    cellregion!(b, 1)
    regionpoint!(b, 5, 5)

    cellregion!(b, 2)
    regionpoint!(b, 5, -1)

    holepoint!(b, 5, 2)

    bregions!(b, grid1, 1:6)
    grid2 = simplexgrid(b; maxvolume = 0.6)
    grid2 = glue(grid1, grid2)
end
#
# ![](glue_2d.png)
#



# Plot generation
using GridVisualize
function generateplots(picdir; Plotter = nothing)
    if isdefined(Plotter, :Makie)
        size = (600, 300)
        Plotter.activate!(; type = "png", visible = false)

        p = builderplot(triangulation_of_domain(); Plotter, size)
        Plotter.save(joinpath(picdir, "triangulation_of_domain.png"),p)
        
        p = builderplot(nicer_triangulation_of_domain(); Plotter, size)
        Plotter.save(joinpath(picdir, "nicer_triangulation_of_domain.png"),p)
        
        p = builderplot(triangulation_of_domain_with_subregions(); Plotter, size)
        Plotter.save(joinpath(picdir, "triangulation_of_domain_with_subregions.png"),p)

        p = builderplot(square_localref(); Plotter, size)
        Plotter.save(joinpath(picdir, "square_localref.png"),p)
        
        p = gridplot(direct_square(); Plotter, size)
        Plotter.save(joinpath(picdir, "direct_square.png"),p)
        
        p = builderplot(swiss_cheese_2d(); Plotter, size)
        Plotter.save(joinpath(picdir, "swiss_cheese_2d.png"),p)
        
        p = gridplot(remesh_2d(); Plotter, size)
        Plotter.save(joinpath(picdir, "remesh_2d.png"),p)

        p = gridplot(glue_2d(); Plotter, size)
        Plotter.save(joinpath(picdir, "glue_2d.png"),p)
    end
end
