"""
$(TYPEDSIGNATURES)

Create Grid from Triangle input data.

See the documentations for 
[`TriangulateIO`](https://juliageometry.github.io/Triangulate.jl/stable/#Triangulate.TriangulateIO),
[`triunsuitable`](https://juliageometry.github.io/Triangulate.jl/stable/#Triangulate.triunsuitable-Tuple{Function})
and the [short](https://juliageometry.github.io/Triangulate.jl/stable/#Triangulate.triangulate-Tuple{String,TriangulateIO})
resp. [long](https://juliageometry.github.io/Triangulate.jl/stable/triangle-h/)  documentation of the Triangle
control flags.
"""
function ExtendableGrids.simplexgrid(flags::String, input::Triangulate.TriangulateIO;unsuitable=nothing)

    if typeof(unsuitable)!=Nothing
        triunsuitable(unsuitable)
    end
        
    triout,vorout=Triangulate.triangulate(flags,input)

    pointlist=triout.pointlist
    if eltype(pointlist)!=Float64
        pointlist=Array{Float64,2}(pointlist)
    end
    
    trianglelist=triout.trianglelist
    if eltype(trianglelist)!=Int32
        trianglelist=Array{Int32,2}(trianglelist)
    end

    cellregions=Vector{Int32}(vec(triout.triangleattributelist))
    
    segmentlist=triout.segmentlist
    if eltype(segmentlist)!=Int32
        segmentlist=Array{Int32,2}(segmentlist)
    end
    
    segmentmarkerlist=triout.segmentmarkerlist
    if eltype(segmentmarkerlist)!=Int32
        segmentmarkerlist=Array{Int32,2}(segmentmarkerlist)
    end

    ExtendableGrids.simplexgrid(pointlist,trianglelist,cellregions,segmentlist,segmentmarkerlist)
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
function triangulateio(;flags::String="pAaqDQ",
                       points=Array{Cdouble,2}(undef,0,0),
                       bfaces=Array{Cint,2}(undef,0,0),
                       bfaceregions=Array{Cint,1}(undef,0),
                       regionpoints=Array{Cdouble,2}(undef,0,0),
                       regionnumbers=Array{Cint,1}(undef,0),
                       regionvolumes=Array{Cdouble,1}(undef,0)
                       )
    @assert ndims(points)==2
    if size(points,2)==2
        points=transpose(points)
    end
    if typeof(points)!=Array{Cdouble,2}
        points=Array{Cdouble,2}(points)
    end
    @assert(size(points,2)>2)
    
    @assert ndims(bfaces)==2
    if size(bfaces,2)==2
        bfaces=transpose(bfaces)
    end
    if typeof(bfaces)!=Array{Cint,2}
        bfaces=Array{Cint,2}(bfaces)
    end
    # @assert(size(bfaces,2)>0)
    
    @assert ndims(bfaceregions)==1
    @assert size(bfaceregions,1)==size(bfaces,2)
    if typeof(bfaceregions)!=Array{Cint,1}
        bfaceregions=Array{Cint,1}(bfaceregions)
    end
    
    @assert ndims(regionpoints)==2
    if size(regionpoints,1)!=2
        regionpoints=transpose(regionpoints)
    end
    if typeof(regionpoints)!=Array{Cdouble,2}
        regionpoints=Array{Cdouble,2}(regionpoints)
    end
    # @assert(size(regionpoints,2)>0)
    
    @assert ndims(regionnumbers)==1
    @assert ndims(regionvolumes)==1
    @assert size(regionnumbers,1)==size(regionpoints,2)
    @assert size(regionvolumes,1)==size(regionpoints,2)
    
    nholes=0
    nregions=0
    for i=1:length(regionnumbers)
        if regionnumbers[i]==0
            nholes+=1
        else
            nregions+=1
        end
    end


    
    regionlist=Array{Cdouble,2}(undef,4,nregions)
    holelist=Array{Cdouble,2}(undef,2,nholes)
    
    ihole=1
    iregion=1
    for i=1:length(regionnumbers)
        if regionnumbers[i]==0
            holelist[1,ihole]=regionpoints[1,i]
            holelist[2,ihole]=regionpoints[2,i]
            ihole+=1
        else
            regionlist[1,iregion]=regionpoints[1,i]
            regionlist[2,iregion]=regionpoints[2,i]
            regionlist[3,iregion]=regionnumbers[i]
            regionlist[4,iregion]=regionvolumes[i]
            iregion+=1
        end
    end
    tio=Triangulate.TriangulateIO()
    tio.pointlist=points
    if size(bfaces,2)>0
        tio.segmentlist=bfaces
    end
    if size(bfaceregions,1)>0
        tio.segmentmarkerlist=bfaceregions
    end
    if size(regionlist,2)>0
        tio.regionlist=regionlist
    end
    if size(holelist,2)>0
        tio.holelist=holelist
    end
    tio
end