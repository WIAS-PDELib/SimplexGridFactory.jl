"""
$(TYPEDSIGNATURES)

Create Grid from TetGen data.

See [`default_options`](@ref) for available `kwargs`.

"""
function ExtendableGrids.simplexgrid(::Type{TetGenType}, TetGen, input; kwargs...)
    opts = blendoptions!(default_options(:tetgen); kwargs...)

    flags = makeflags(opts, :tetgen)

    if opts[:verbose]
        @show flags
    end

    if !isnothing(opts[:unsuitable])
        TetGen.tetunsuitable!(opts[:unsuitable])
    end

    tetout = TetGen.tetrahedralize(input, flags)

    return ExtendableGrids.simplexgrid(tetout)
end

"""
$(TYPEDSIGNATURES)
Create a RawTetGenIO structure 
from a number of input arrays.
The 2D input arrays are transposed if necessary and converted to
the proper data types for TetGen.
 
This conversion is not performed if the data types are those
indicated in the defaults and the leading dimension of 2D arrays
corresponds to the space dimension.

"""
function tetgenio(
        TetGen; points = Array{Cdouble, 2}(undef, 0, 0),
        bfaces = Array{Cint, 2}(undef, 0, 0),
        bfaceregions = Array{Cint, 1}(undef, 0),
        regionpoints = Array{Cdouble, 2}(undef, 0, 0),
        regionnumbers = Array{Cint, 1}(undef, 0),
        regionvolumes = Array{Cdouble, 1}(undef, 0)
    )
    @assert ndims(points) == 2
    if size(points, 2) == 3
        points = transpose(points)
    end
    if typeof(points) != Array{Cdouble, 2}
        points = Array{Cdouble, 2}(points)
    end
    @assert(size(points, 2) > 2)

    # if  ndims(bfaces)==2
    #     if size(bfaces,2)==2
    #         bfaces=transpose(bfaces)
    #     end
    #     if typeof(bfaces)!=Array{Cint,2}
    #         bfaces=Array{Cint,2}(bfaces)
    #     end
    # end
    # @assert(size(bfaces,2)>0)

    @assert ndims(bfaceregions) == 1
    if ndims(bfaces) == 2
        @assert size(bfaceregions, 1) == size(bfaces, 2)
    else
        @assert size(bfaceregions, 1) == size(bfaces, 1)
    end
    if typeof(bfaceregions) != Array{Cint, 1}
        bfaceregions = Array{Cint, 1}(bfaceregions)
    end

    @assert ndims(regionpoints) == 2
    if size(regionpoints, 1) != 3
        regionpoints = transpose(regionpoints)
    end
    if typeof(regionpoints) != Array{Cdouble, 2}
        regionpoints = Array{Cdouble, 2}(regionpoints)
    end
    # @assert(size(regionpoints,2)>0)

    @assert ndims(regionnumbers) == 1
    @assert ndims(regionvolumes) == 1
    @assert size(regionnumbers, 1) == size(regionpoints, 2)
    @assert size(regionvolumes, 1) == size(regionpoints, 2)

    nholes = 0
    nregions = 0
    for i in 1:length(regionnumbers)
        if regionnumbers[i] == 0
            nholes += 1
        else
            nregions += 1
        end
    end

    regionlist = Array{Cdouble, 2}(undef, 5, nregions)
    holelist = Array{Cdouble, 2}(undef, 3, nholes)

    ihole = 1
    iregion = 1
    for i in 1:length(regionnumbers)
        if regionnumbers[i] == 0
            holelist[1, ihole] = regionpoints[1, i]
            holelist[2, ihole] = regionpoints[2, i]
            holelist[3, ihole] = regionpoints[3, i]
            ihole += 1
        else
            regionlist[1, iregion] = regionpoints[1, i]
            regionlist[2, iregion] = regionpoints[2, i]
            regionlist[3, iregion] = regionpoints[3, i]
            regionlist[4, iregion] = regionnumbers[i]
            regionlist[5, iregion] = regionvolumes[i]
            iregion += 1
        end
    end
    tio = TetGen.RawTetGenIO{Float64}()
    tio.pointlist = points
    if size(bfaces, 2) > 0
        TetGen.facetlist!(tio, bfaces)
    end
    if size(bfaceregions, 1) > 0
        tio.facetmarkerlist = bfaceregions
    end
    if size(regionlist, 2) > 0
        tio.regionlist = regionlist
    end
    if size(holelist, 2) > 0
        tio.holelist = holelist
    end
    return tio
end

"""
$(TYPEDSIGNATURES)

Create tetgen input from the current state of the builder.
"""
function tetgenio(this::SimplexGridBuilder)
    dim_space(this) = 3 || throw(error("dimension !=2 not implemented"))

    return tetgenio(
        this.Generator; points = this.pointlist.points,
        bfaces = this.facets,
        bfaceregions = this.facetregions,
        regionpoints = this.regionpoints,
        regionnumbers = this.regionnumbers,
        regionvolumes = this.regionvolumes
    )
end

function writesmesh(io::IO, tio)
    shift = 1
    npoints = size(tio.pointlist, 2)
    dimension = 3
    npointboundarymarkers = 0
    npointattributes = 0
    nfacets = length(tio.facetlist)
    nfacetboundarymarkers = length(tio.facetmarkerlist) > 0 ? 1 : 0
    nholes = size(tio.holelist, 2)
    nregions = size(tio.regionlist, 2)

    write(io, @sprintf("# TetGen input file\n"))
    write(io, @sprintf("# Created by SimplexGridFactory.jl %s\n", string(now())))
    write(io, @sprintf("\n# part 1: node list."))
    write(io, @sprintf("\n%ld %d %d %d", npoints, dimension, npointboundarymarkers, npointattributes))
    for ipoint in 1:npoints
        write(
            io, @sprintf(
                "\n%ld %.17g %.17g %.17g", ipoint - 1,
                tio.pointlist[1, ipoint], tio.pointlist[2, ipoint], tio.pointlist[3, ipoint]
            )
        )
    end
    write(io, @sprintf("\n\n# part 2: facet list."))
    write(io, @sprintf("\n%ld %d", nfacets, nfacetboundarymarkers))
    for ifacet in 1:nfacets
        facet = tio.facetlist[ifacet].polygonlist[1]
        ncorners = length(facet)
        write(io, @sprintf("\n%ld ", ncorners))
        for icorner in 1:ncorners
            write(io, @sprintf("%ld ", facet[icorner] - 1))
        end
        if nfacetboundarymarkers > 0
            write(io, @sprintf("%ld ", tio.facetmarkerlist[ifacet]))
        end
    end
    write(io, @sprintf("\n\n# part 3: hole list."))
    write(io, @sprintf("\n%ld", nholes))
    for ihole in 1:nholes
        write(
            io, @sprintf(
                "\n%ld %.17g %.17g %.17g", ihole - 1,
                tio.holelist[1, ihole], tio.holelist[2, ihole], tio.holelist[3, ihole]
            )
        )
    end
    write(io, @sprintf("\n\n# part 4: region list."))
    write(io, @sprintf("\n%ld", nregions))
    for iregion in 1:nregions
        write(
            io, @sprintf(
                "\n%ld %.17g %.17g %.17g %f %.17g", iregion - 1,
                tio.regionlist[1, iregion], tio.regionlist[2, iregion], tio.regionlist[3, iregion],
                tio.regionlist[4, iregion], tio.regionlist[5, iregion]
            )
        )
    end
    write(io, "\n")
    return
end

function writesmesh(fname::String, tetin)
    open(fname, "w") do io
        writesmesh(io, tetin)
    end
    return
end
