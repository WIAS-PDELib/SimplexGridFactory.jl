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

        center = zeros(2)

        new_unsuitable = (x1, y1, x2, y2, x3, y3, area) -> begin

            # center of the triangle
            center[1] = (x1 + x2 + x3) / 3
            center[2] = (y1 + y2 + y3) / 3

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
        center = zeros(3)
        matrix = @MArray zeros(3, 3)

        new_unsuitable = (pa, pb, pc, pd) -> begin
            @. center = 0.25 * (pa + pb + pc + pd)
            @. matrix[:, 1] = pb - pa
            @. matrix[:, 2] = pc - pa
            @. matrix[:, 3] = pd - pa

            vol = abs(det(matrix)) / 6
            return unsuitable(center, vol)
        end
    else
        new_unsuitable = unsuitable
    end

    return Generator.tetunsuitable!(new_unsuitable)

end
