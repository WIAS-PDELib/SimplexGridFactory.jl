"""
    istetgen(Generator)

Check if generator is TetGen.
"""
istetgen(g) = typeof(g) == Module && nameof(g) == :TetGen

"""
    istriangulate(Generator)

Check if generator is Triangulate.
"""
istriangulate(g) = typeof(g) == Module && nameof(g) == :Triangulate

abstract type TetGenType end

abstract type TriangulateType end

function set_unsuitable!(::Type{TriangulateType}, Generator, unsuitable)

    # We assume that the "Triangle" uses Cdouble internally
    if hasmethod(unsuitable, Tuple{Vector{Float64}, Float64})
        new_unsuitable = (x1, y1, x2, y2, x3, y3, area) -> begin
            center = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3] # center of the triangle
            return unsuitable(center, area)
        end
    else
        new_unsuitable = unsuitable
    end

    return Generator.triunsuitable!(new_unsuitable)

end

function set_unsuitable!(::Type{TetGenType}, Generator, unsuitable)

    # We assume that the "TetGen" uses Cdouble internally
    if hasmethod(unsuitable, Tuple{Vector{Float64}, Float64})

        new_unsuitable = (pa, pb, pc, pd) -> begin
            center = 0.25 * (pa + pb + pc + pd)
            vol = abs(det(hcat(pb - pa, pc - pa, pd - pa))) / 6
            return unsuitable(center, vol)
        end
    else
        new_unsuitable = unsuitable
    end

    return Generator.tetunsuitable!(new_unsuitable)

end
