"""
$(TYPEDEF)
    
Simplex grid builder: wrapper around array based mesh generator interface.
It allows to build up the input data incrementally.

"""
mutable struct SimplexGridBuilder
    Generator::Union{Module, Nothing}
    current_facetregion::Cint
    current_cellregion::Cint
    current_cellvolume::Cdouble
    point_identity_tolerance::Cdouble
    facetregions::Vector{Cint}
    facets::Vector{Vector{Cint}}
    pointlist::BinnedPointList{Cdouble}
    regionpoints::ElasticArray{Cdouble, 2}
    regionnumbers::Vector{Cint}
    regionvolumes::Vector{Cdouble}
    options::Dict{Symbol, Any}
    checkexisting::Bool
    _savedpoint::Cint
    SimplexGridBuilder(x::Nothing) = new()
end

"""
```
SimplexGridBuilder(; Generator=nothing,
                     tol=1.0e-12,
                     checkexisting=true)
```

Create a SimplexGridBuilder.

- `Generator`: module corresponding to mesh generator package.
   Valid choices are `TetGen` and `Triangulate`, corresponding to the 
   respective Julia packages.
- `checkexisting`: whether to check for already existing points
- `tol`: two points below this tolerance will be merged if `checkexisting` is true
"""
function SimplexGridBuilder(; Generator = nothing, tol = 1.0e-12, checkexisting = true)
    builder = SimplexGridBuilder(nothing)

    if isnothing(Generator)
        throw(ArgumentError("Missing Generator: SimplexGridBuilder needs Generator=TetGen or Generator=Triangulate as argument"))
    end

    builder.Generator = Generator

    if istetgen(Generator)
        dim = 3
        mesher = :tetgen
    elseif istriangulate(Generator)
        dim = 2
        mesher = :triangle
    else
        throw(ArgumentError("Wrong Generator: SimplexGridBuilder needs Generator=TetGen or Generator=Triangulate as argument"))
    end

    builder.current_facetregion = 1
    builder.current_cellregion = 1
    builder.current_cellvolume = 1

    builder.point_identity_tolerance = tol
    builder.facets = []
    builder.facetregions = []
    builder.pointlist = BinnedPointList(dim; tol = tol)
    builder.regionpoints = ElasticArray{Cdouble}(undef, dim, 0)
    builder.regionvolumes = []
    builder.regionnumbers = []
    builder.options = default_options(mesher)
    builder.checkexisting = checkexisting
    builder._savedpoint = 0
    return builder
end

"""
$(SIGNATURES)

Whether to check for already existing points
"""
checkexisting!(builder, b) = builder.checkexisting = b

"""
$(SIGNATURES)

Set some mesh generation options, see [`default_options`](@ref)
"""
options!(builder::SimplexGridBuilder; kwargs...) = blendoptions!(builder.options; kwargs...)

"""
$(SIGNATURES)

Space dimension of builder.
"""
ExtendableGrids.dim_space(builder::SimplexGridBuilder) = size(builder.pointlist.points, 1)

"""
```
point!(builder,x)
point!(builder,x,y)
point!(builder,x,y,z)
point!(builder,vec_or_tuple)
```

Add point or merge with already existing point. Returns its index
which can be used to set up facets with [`facet!`](@ref).
"""
function point!(builder::SimplexGridBuilder, x)
    dim_space(builder) == 1 || throw(DimensionMismatch())
    return insert!(builder.pointlist, [x])
end

function point!(builder::SimplexGridBuilder, x, y)
    dim_space(builder) == 2 || throw(DimensionMismatch())
    return insert!(builder.pointlist, [x, y])
end

function point!(builder::SimplexGridBuilder, x, y, z)
    dim_space(builder) == 3 || throw(DimensionMismatch())
    return insert!(builder.pointlist, [x, y, z])
end

const PointCoord = Union{AbstractVector, Tuple}

point!(builder::SimplexGridBuilder, p::PointCoord) = point!(builder, p...)

"""
```
cellregion!(builder,region)
```

Set the current cell region (acts on subsequent regionpoint() calls)

Cell regions can be used to distinguish cells of different materials etc.
In the API they are characterized by
- region number set via `cellregion!`
- maximum cell volume set via [`maxvolume!`](@ref)
- region point set via  [`regionpoint!`](@ref). This is some point located
  within the respective region which must be surrounded by facets in a watertight 
  manner. 
"""
cellregion!(builder::SimplexGridBuilder, i) = builder.current_cellregion = i

"""
```
maxvolume!(builder,vol)
```

Set the current cell volume resp. area (acts on subsequent regionpoint() calls).
See [`cellregion!`](@ref).
"""
maxvolume!(builder::SimplexGridBuilder, vol) = builder.current_cellvolume = vol

"""
```
regionpoint!(builder,x)
regionpoint!(builder,x,y)
regionpoint!(builder,x,y,z)
regionpoint!(builder,vec_or_tuple)
```
Add a region point marking a region, using current cell volume an cell region
See [`cellregion!`](@ref).
"""
function regionpoint!(builder::SimplexGridBuilder, x)
    dim_space(builder) == 1 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x))
    push!(builder.regionvolumes, builder.current_cellvolume)
    return push!(builder.regionnumbers, builder.current_cellregion)
end

function regionpoint!(builder::SimplexGridBuilder, x, y)
    dim_space(builder) == 2 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x, y))
    push!(builder.regionvolumes, builder.current_cellvolume)
    return push!(builder.regionnumbers, builder.current_cellregion)
end

function regionpoint!(builder::SimplexGridBuilder, x, y, z)
    dim_space(builder) == 3 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x, y, z))
    push!(builder.regionvolumes, builder.current_cellvolume)
    return push!(builder.regionnumbers, builder.current_cellregion)
end

regionpoint!(builder::SimplexGridBuilder, p::PointCoord) = regionpoint!(builder, p...)

"""
```
holepoint!(builder,x)
holepoint!(builder,x,y)
holepoint!(builder,x,y,z)
holepoint!(builder,vec_or_tuple)
```
Add a point marking a hole region. Hole regions need to be surrounded by facets
in a watertight manner.
"""
function holepoint!(builder::SimplexGridBuilder, x)
    dim_space(builder) == 1 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x))
    push!(builder.regionvolumes, 0)
    push!(builder.regionnumbers, 0)
    return nothing
end

function holepoint!(builder::SimplexGridBuilder, x, y)
    dim_space(builder) == 2 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x, y))
    push!(builder.regionvolumes, 0)
    push!(builder.regionnumbers, 0)
    return nothing
end

function holepoint!(builder::SimplexGridBuilder, x, y, z)
    dim_space(builder) == 3 || throw(DimensionMismatch())
    append!(builder.regionpoints, (x, y, z))
    push!(builder.regionvolumes, 0)
    push!(builder.regionnumbers, 0)
    return nothing
end

holepoint!(builder::SimplexGridBuilder, p::PointCoord) = holepoint!(builder, p...)

"""
```
facetregion!(builder,region)
```
Set the current facet region. Subsequent facets will be marked with this number.
Facet regions can be used to mark different parts of the boundary, e.g. for
distinguishing boundary conditions.
"""
facetregion!(builder::SimplexGridBuilder, i) = builder.current_facetregion = i

"""
```
facet!(builder,i1)
facet!(builder,i1,i2)
facet!(builder,i1,i2,i3,i4)
facet!(builder,vector_or_tuple)
facet!(builder, (x1,y1), (x2,y2))
facet!(builder, (x1,y1,z1), (x2,y2,z2),(x3,y3,z3))
```

Add a facet via the corresponding point indices returned
by [`point!`](@ref). 

Facets of two points are solely used for 2D grids. Facets
with more than two points are used for 3D grids and must be 
planar.
"""
function facet!(builder::SimplexGridBuilder, i)
    dim_space(builder) == 1 || throw(DimensionMismatch())
    push!(builder.facets, [i])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

function facet!(builder::SimplexGridBuilder, i1, i2)
    dim_space(builder) == 2 || throw(DimensionMismatch())
    push!(builder.facets, [i1, i2])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

function facet!(builder::SimplexGridBuilder, i1, i2, i3)
    dim_space(builder) == 3 || throw(DimensionMismatch())
    push!(builder.facets, [i1, i2, i3])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

function facet!(builder::SimplexGridBuilder, i1, i2, i3, i4)
    dim_space(builder) == 3 || throw(DimensionMismatch())
    push!(builder.facets, [i1, i2, i3, i4])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

function facet!(builder::SimplexGridBuilder, p::Union{Vector, Tuple})
    if dim_space(builder) == 1
        length(p) == 1 || throw(DimensionMismatch())
    end
    if dim_space(builder) == 2
        length(p) == 2 || throw(DimensionMismatch())
    end
    if dim_space(builder) == 3
        length(p) >= 3 || throw(DimensionMismatch())
    end
    push!(builder.facets, [p...])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

facet!(builder::SimplexGridBuilder, p1::PointCoord, p2::PointCoord) = facet!(builder, point!(builder, p1), point!(builder, p2))
function facet!(builder::SimplexGridBuilder, p1::PointCoord, p2::PointCoord, p3::PointCoord)
    return facet!(builder, point!(builder, p1), point!(builder, p2), point!(builder, p3))
end

"""
```
polyfacet!(builder,vector_or_tuple)
```

Add a polygonal facet via the corresponding point indices returned
by [`point!`](@ref). 

Facets with more than two points are used for 3D grids and must be 
planar.
"""
function polyfacet!(builder::SimplexGridBuilder, p::Union{Vector, Tuple})
    push!(builder.facets, [p...])
    push!(builder.facetregions, builder.current_facetregion)
    return length(builder.facets)
end

"""
``` 
simplexgrid(builder; kwargs...)
```

Build simplex grid from the current state of the builder.
`kwargs` overwrite those set with the [`options!`](@ref) method.
See [`default_options`](@ref) for available `kwargs`.
"""
function ExtendableGrids.simplexgrid(builder::SimplexGridBuilder; kwargs...)
    if dim_space(builder) == 2
        facets = Array{Cint, 2}(undef, 2, length(builder.facets))
        for i in 1:length(builder.facets)
            facets[1, i] = builder.facets[i][1]
            facets[2, i] = builder.facets[i][2]
        end
        make_tio = triangulateio
        generator_type = TriangulateType
    else
        facets = builder.facets
        make_tio = tetgenio
        generator_type = TetGenType
    end

    options = blendoptions!(copy(builder.options); kwargs...)

    tio = make_tio(
        builder.Generator;
        points = builder.pointlist.points,
        bfaces = facets,
        bfaceregions = builder.facetregions,
        regionpoints = builder.regionpoints,
        regionnumbers = builder.regionnumbers,
        regionvolumes = builder.regionvolumes
    )

    return ExtendableGrids.simplexgrid(generator_type, builder.Generator, tio; options...)
end

"""
    flags(builder)

Return mesh generator specific flag string created from builder options.
"""
function flags(builder::SimplexGridBuilder)
    return if istetgen(builder.Generator)
        makeflags(builder.options, :tetgen)
    elseif istriangulate(builder.Generator)
        makeflags(builder.options, :triangle)
    else
        nothing
    end
end

"""
     maybewatertight(this::SimplexGridBuilder; bregions=nothing)

Check if facets belonging to boundare regions in bregions are watertight.
This is based on a number of heuristics, only a negative
answer is definitive.
"""
function maybewatertight(this::SimplexGridBuilder; bregions = nothing)
    maybe = true

    bfaceregions = this.facetregions
    if bregions == nothing
        bregions = unique(bfaceregions)
    end

    @info "Checking for dangling facets"
    points = this.pointlist.points
    facets = this.facets
    dim = size(points, 1)
    npoints = size(points, 2)
    nfacets = size(facets, 1)
    ptmarkers = zeros(Int, npoints)

    for ifacet in 1:nfacets
        if bfaceregions[ifacet] ∈ bregions
            for idim in 1:dim
                ptmarkers[facets[ifacet][idim]] += 1
            end
        end
    end

    uptmarkers = unique(ptmarkers)
    if 1 ∈ uptmarkers
        maybe = false
    end

    if dim == 3 && 2 ∈ uptmarkers
        maybe = false
    end

    if maybe
        @info "Maybe description is watertight, but not sure"
    else
        @warn "Description is not watertight"
        for ifacet in 1:nfacets
            if bfaceregions[ifacet] ∈ bregions
                for idim in 1:dim
                    pt = facets[ifacet][idim]
                    if ptmarkers[pt] < dim
                        @warn "Dangling facet $ifacet (bregion $(bfaceregions[ifacet]), point $(points[:, pt])"
                    end
                end
            end
        end
    end
    return maybe
end
