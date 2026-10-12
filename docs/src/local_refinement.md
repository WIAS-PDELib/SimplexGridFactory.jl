## Local refinement

### 2D - Triangulate

The `unsuitable` mechanism is used to refine the 2D triangulation on a domain, for instance surrounding a geometric object (e.g., a point or line). 
Refinement is triggered when a condition for a given triangle is satisfied.
For expamle when the distance between a point of interest and the center (barycentric coordinates) of the triangle falls below a threshold defined by the triangle’s area.
#### Example usage:

- One could define the `unsuitable` function with 6 points and the area of the triangle (x1, y1, x2, y2, x3, y3, area). This is the default signature in [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl):
```julia
    refinement_center = [0.5, 0.5]
    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y2 + y2 + y3) / 3]
        dist = norm(bary - refinement_center)
        return area > 0.01 * dist ? 1 : 0
    end
```  

- As a simpler alternative, one could directly use the barycentric coordinates and the area of the triangle (center, area): 
```julia
    refinement_center = [0.5, 0.5]
    function unsuitable(center, area)
        dist = norm(center - refinement_center)
        return area > 0.01 * dist ? 1 : 0
    end
```
To hand this function over to the builder:     
```julia 
    options!(builder; unsuitable = unsuitable)
```

Full example can be found [here](../script_examples/examples2d/#Local-refinement)

### 3D - TetGen

The  `unsuitable` mechanism is also used to refine the 3D triangulation on a domain, for instance surrounding a geometric object (e.g., a point or line). Refinement is triggered when a condition for a given tetrahedron is satisfied.
For example when the distance between the point of interest and the barycentric coordinates of a tetrahedron falls below a threshold defined by the tetrahedron’s volume.   

#### Example usage:

- One could define the `unsuitable` function with 4 points (pa, pb, pc, pd). 
This is the default signature in [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl):  
```julia
    refinement_center = [0.5, 0.5, 0.5]
    function unsuitable(pa, pb, pc, pd)
        vol = det(hcat(pb - pa, pc - pa, pd - pa)) / 6
        center = 0.25 * (pa + pb + pc + pd) 
        return vol > 0.05 * norm(center - refinement_center)^2.5
    end
```


- As a simpler alternative, one could directly use the barycentric coordinates and the volume of the tetrahedron (center, volume):  
```julia
    refinement_center = [0.5, 0.5, 0.5]
    function unsuitable(center, vol)
        return vol > 0.05 * norm(center - refinement_center)^2.5
    end

```
   To hand this function over to the builder:     
```julia 
    options!(builder; unsuitable = unsuitable)
```
