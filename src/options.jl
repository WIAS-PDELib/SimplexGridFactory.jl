"""
$(SIGNATURES)
Create dictionary of mesh generation options with default values. These at once describe the
keyword arguments available to the methods of the package and are listed in the following table:


| keyword        | default | 2D  | 3D  | Explanation                                                |
|:---------------|:-------:|:---:|:---:|:-----------------------------------------------------------|
| PLC            | true    | -p  | -p  | Triangulate/tetraheralize PLSG/PLC                         |
| refine         | false   | -r  | -r  | Refines a previously generated mesh.                       |
| quality        | true    | -q  | -q  | Quality mesh generation                                    |
| minangle       | 2D: 20  |     |     | Minimum angle for quality                                  |
|                | 3D: 10  |     |     |                                                            |
| volumecontrol  | true    | -a  | -a  | Maximum area constraint                                    |
| maxvolume      | Inf     |     |     | Value of area/volume constraint if less than Inf           |
| attributes     | true    | -A  | -A  | Regional attribute to each simplex.                        |
| confdelaunay   | true    | -D  |     | Ensure that all circumcenter lie within the domain.        |
| nosteiner      | false   | -Y  | -Y  | Prohibits insertion of Steiner points on the mesh boundary |
| quiet          | true    | -Q  | -Q  | Suppress all output unless an error occurs.                |
| verbose        | false   | -V  | -V  | Give detailed information.                                 |
| debugfacets    | true    |     | -d  | Detects self-intersections of facets of the PLC.           |
| check          | false   | -C  | -C  | Checks the consistency of the final mesh.                  |
| optlevel       | 1       |     | -O  | Specifies the level of mesh optimization.                  |
| unsuitable     | nothing |     |     | Unsuitable function                                        |
| addflags       | ""      |     |     | Additional flags                                           |
| flags          | nothing |     |     | Set flags, overwrite all other options                     |


For mesh generation, these are turned into mesh generator control flags. This process can be completely
overwritten by specifying the flags parameter. 

For the flags parameter in 2D see the
[short](https://juliageometry.github.io/Triangulate.jl/stable/#Triangulate.triangulate-Tuple{String,TriangulateIO})
resp. [long](https://juliageometry.github.io/Triangulate.jl/stable/triangle-h/)
documentation of the Triangle control flags.

For the 3D case, see the corresponding [TetGen flags](https://juliageometry.github.io/TetGen.jl/stable/#TetGen.tetrahedralize-Tuple{RawTetGenIO{Float64},String})

The `unsuitable` parameter should be a function, see
[`triunsuitable!`](https://juliageometry.github.io/TetGen.jl/stable/#TetGen.triunsuitable!-Tuple{Function}) .

"""
default_options(mesher) = Dict{Symbol, Any}(
    :PLC => true,
    :refine => false,
    :quality => true,
    :minangle => mesher == :triangle ? 20 : 10,
    :volumecontrol => true,
    :maxvolume => Inf,
    :attributes => true,
    :confdelaunay => true,
    :optlevel => 1,
    :nosteiner => false,
    :quiet => true,
    :verbose => false,
    :debugfacets => true,
    :check => false,
    :unsuitable => nothing,
    :flags => nothing,
    :addflags => ""
)

function blendoptions!(opt; kwargs...)
    for (k, v) in kwargs
        if haskey(opt, Symbol(k))
            opt[Symbol(k)] = v
        end
    end

    if opt[:verbose]
        for (k, v) in kwargs
            if !haskey(opt, Symbol(k))
                println("Warning: ignored kwarg $(k)=$(v) for simplexgrid")
            end
        end
    end
    return opt
end

function makeflags(options, mesher)
    return if isnothing(options[:flags])
        flags = ""
        options[:PLC] ? flags *= "p" : nothing
        options[:refine] ? flags *= "r" : nothing
        if options[:quality]
            minangle = Float64(options[:minangle])
            if mesher == :tetgen
                flags *= @sprintf("q/%.2f", minangle)
            else
                flags *= @sprintf("q%.2f", minangle)
            end
        end
        if options[:volumecontrol]
            flags *= "a"
            maxvolume = options[:maxvolume]
            if !isinf(maxvolume)
                flags *= @sprintf("%.40f", maxvolume)
            end
        end
        options[:attributes] ? flags *= "A" : nothing
        !isnothing(options[:unsuitable]) && mesher == :triangle ? flags *= "u" : nothing
        options[:confdelaunay] ? flags *= "D" : nothing
        if options[:optlevel] > 0 && mesher == :tetgen
            optlevel = options[:optlevel]
            flags *= @sprintf("O%d", optlevel)
        end
        options[:nosteiner] ? flags *= "Y" : nothing
        options[:quiet] ? flags *= "Q" : nothing
        options[:verbose] ? flags *= "V" : nothing
        options[:debugfacets] && mesher == :teten ? flags *= "d" : nothing
        options[:check] ? flags *= "C" : nothing
        flags *= options[:addflags]
    else
        options[:flags]
    end
end
