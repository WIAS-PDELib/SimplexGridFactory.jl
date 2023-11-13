"""
$(TYPEDSIGNATURES)

Create Grid from Triangle input data.

See [`default_options`](@ref) for available `kwargs`.
"""
function ExtendableGrids.simplexgrid(::Type{TriangulateType}, Triangulate, input; kwargs...)
    opts = blendoptions!(default_options(); kwargs...)

    flags = makeflags(opts, :triangle)

    if opts[:verbose]
        @show flags
    end
    if !isnothing(opts[:unsuitable])
        Triangulate.triunsuitable(opts[:unsuitable])
    end

    triout, vorout = Triangulate.triangulate(flags, input)

    pointlist = triout.pointlist

    trianglelist = triout.trianglelist

    if size(triout.triangleattributelist, 2) == 0
        # Add default for cellregions if that was not created
        cellregions = ones(Int32, size(trianglelist, 2))
    else
        cellregions = Vector{Int32}(vec(triout.triangleattributelist))
    end

    segmentlist = triout.segmentlist

    segmentmarkerlist = triout.segmentmarkerlist

    if size(pointlist, 2) == 0
        error("Emtpy list of generated points. May be the geometry description is not watertight ?") |> throw
    end

    if size(trianglelist, 2) == 0
        error("Emtpy list of generated triangles. May be the geometry description is not watertight ?") |> throw
    end

    ExtendableGrids.simplexgrid(pointlist, trianglelist, cellregions, segmentlist, segmentmarkerlist)
end

"""
$(TYPEDSIGNATURES)
Create a TriangulateIO structure 
from a number of input arrays.
The 2D input arrays are transposed if necessary and converted to
the proper data types for Triangulate.
 
This conversion is not performed if the data types are those
indicated in the defaults and the leading dimension of 2D arrays
corresponds to the space dimension.

"""
function triangulateio(Triangulate;
                       points = Array{Cdouble, 2}(undef, 0, 0),
                       bfaces = Array{Cint, 2}(undef, 0, 0),
                       bfaceregions = Array{Cint, 1}(undef, 0),
                       regionpoints = Array{Cdouble, 2}(undef, 0, 0),
                       regionnumbers = Array{Cint, 1}(undef, 0),
                       regionvolumes = Array{Cdouble, 1}(undef, 0))
    ndims(points) == 2 || throw(DimensionMismatch("Wrong space dimension"))
    if size(points, 2) == 2
        points = transpose(points)
    end
    if typeof(points) != Array{Cdouble, 2}
        points = Array{Cdouble, 2}(points)
    end
    size(points, 2) > 2 || throw(ErrorException("Expected more than 2 input points"))

    @assert ndims(bfaces) == 2
    if size(bfaces, 2) == 2
        bfaces = transpose(bfaces)
    end
    if typeof(bfaces) != Array{Cint, 2}
        bfaces = Array{Cint, 2}(bfaces)
    end
    # @assert(size(bfaces,2)>0)

    ndims(bfaceregions) == 1 || throw(DimensionMismatch("bfaceregions must be vector"))
    size(bfaceregions, 1) == size(bfaces, 2) || throw(DimensionMismatch("size(bfaceregions,1) != size(bfaces,2)"))
    if typeof(bfaceregions) != Array{Cint, 1}
        bfaceregions = Array{Cint, 1}(bfaceregions)
    end

    ndims(regionpoints) == 2 || throw(DimensionMismatch("Region point must be 2D"))
    if size(regionpoints, 1) != 2
        regionpoints = transpose(regionpoints)
    end
    if typeof(regionpoints) != Array{Cdouble, 2}
        regionpoints = Array{Cdouble, 2}(regionpoints)
    end
    # @assert(size(regionpoints,2)>0)

    @assert ndims(regionnumbers) == 1 || throw(DimensionMismatch("regionnumbers  must be vector"))
    @assert ndims(regionvolumes) == 1 || throw(DimensionMismatch("regionvolumes  must be vector"))
    @assert size(regionnumbers, 1) == size(regionpoints, 2) ||
            throw(DimensionMismatch("size(regionnumbers,1) != size(regionpoints,2)"))
    @assert size(regionvolumes, 1) == size(regionpoints, 2) ||
            throw(DimensionMismatch("size(regionvolumes,1) !== size(regionpoints,2)"))

    nholes = 0
    nregions = 0
    for i = 1:length(regionnumbers)
        if regionnumbers[i] == 0
            nholes += 1
        else
            nregions += 1
        end
    end

    regionlist = Array{Cdouble, 2}(undef, 4, nregions)
    holelist = Array{Cdouble, 2}(undef, 2, nholes)

    ihole = 1
    iregion = 1
    for i = 1:length(regionnumbers)
        if regionnumbers[i] == 0
            holelist[1, ihole] = regionpoints[1, i]
            holelist[2, ihole] = regionpoints[2, i]
            ihole += 1
        else
            regionlist[1, iregion] = regionpoints[1, i]
            regionlist[2, iregion] = regionpoints[2, i]
            regionlist[3, iregion] = regionnumbers[i]
            regionlist[4, iregion] = regionvolumes[i]
            iregion += 1
        end
    end

    tio = Triangulate.TriangulateIO()
    tio.pointlist = points

    if size(bfaces, 2) > 0
        tio.segmentlist = bfaces
    end

    if size(bfaceregions, 1) > 0
        tio.segmentmarkerlist = bfaceregions
    end

    if size(regionlist, 2) > 0
        tio.regionlist = regionlist
    end

    if size(holelist, 2) > 0
        tio.holelist = holelist
    end

    tio
end

"""
$(TYPEDSIGNATURES)

Create triangle input from the current state of the builder.
"""
function triangulateio(this::SimplexGridBuilder)
    dim_space(this) == 2 || throw(error("dimension !=2 not implemented"))
    facets = Array{Cint, 2}(undef, 2, length(this.facets))
    for i = 1:length(this.facets)
        facets[1, i] = this.facets[i][1]
        facets[2, i] = this.facets[i][2]
    end

    triangulateio(this.Generator;
                  points = this.pointlist.points,
                  bfaces = facets,
                  bfaceregions = this.facetregions,
                  regionpoints = this.regionpoints,
                  regionnumbers = this.regionnumbers,
                  regionvolumes = this.regionvolumes)
end
