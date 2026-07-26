### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
    using Revise
    using ExtendableGrids
    using SimplexGridFactory
    using GridVisualize
    using Triangulate
    using PlutoVista
    using VoronoiFVM
    default_plotter!(PlutoVista)
end

# ╔═╡ 732e7d64-0b14-409b-9776-c222951bd79f
md"""
# Creating a grid for a cylinder

We want to have grid refinement towards the wall in order to be eventually able to calculate correct heat transfer through the wall.
"""

# ╔═╡ 84d3f058-3838-4c38-afc5-f2ea30fb9bb6
begin
    nref = 2 # refinement level
    rad = 1 # outer radius
    rad_inner = 0.5 * rad  # inner radius - transition between ring and inner part
    hmin = rad * 0.025 * 2.0^(-nref) # smallest radial grid size (close to wall)
    hmax = rad * 0.25 * 2.0^(-nref) # largest radial grid size
    nang = 15 * 2^(nref) |> ceil |> Int # resolution in angular direction
    len = 10 # length of cylinder
    nz = len * 2^nref + 1 |> ceil |> Int # resolution in length direction
end

# ╔═╡ 08c71069-4ecb-450c-a905-ca6e7a62ac88
md"""
## Step 1: create a ring grid
"""

# ╔═╡ ea4d6975-bbf9-4b72-b3df-aa4d1faea6c7
R = geomspace(rad_inner, rad, hmax, hmin)

# ╔═╡ a5120dd7-4442-4f87-90c0-bdbce949cd31
Φ = range(0, 2π, length = nang)

# ╔═╡ 364db500-712e-4b02-b124-1fa005e53e8e
g_ring = ringsector(R, Φ)

# ╔═╡ f4c5efff-d36a-445c-b59d-abfdc5690212
gridplot(g_ring)

# ╔═╡ 79e6fa42-fdb2-4b99-a3ee-f5281a640fb2
md"""
In order to apply the Voronoi finite volume method we need a boundary conforming
Delaunay grid which ensures that edge transfer coefficients are nonnegative. 
In 2D, this means that for interior edges, the sum of the angels opposite to an edge is less or equal than π. 
For edges situated at an interior or exterior boundary, the angle opposite to the edge should be less or equal than ``\frac{π}{2}``.   

`VoronoiFVM.nondelaunay` returns the edges (as pairs of points) where this coefficient is negative due to a violation of the boundary conforming Delaunay property.
"""

# ╔═╡ 6784f6de-df9d-431b-b4b5-ca23f1bcde65
md"""
The grid created as above misses this condition:
"""

# ╔═╡ 83070579-667c-4475-9716-0980a2447fdf
nondelaunay(g_ring)

# ╔═╡ 436ede6d-7c88-40f1-bf0f-b2691447314f
md"""
In order to fix this situation while ensuring grid refinement towards the 
outer boundary, we re-create the grid using `Triangulate.jl` with the same
point list.
"""

# ╔═╡ 2638116e-d895-4099-a6e4-2a996512cb4a
begin
    b1 = SimplexGridBuilder(Generator = Triangulate)
    coord = g_ring[Coordinates]
    for i in 1:size(coord, 2)
        p = point!(b1, coord[:, i])
    end
    bregions!(b1, g_ring)
    cellregion!(b1, 1)
    regionpoint!(b1, 0, -0.75)
    holepoint!(b1, 0, 0)
    g_ring2 = simplexgrid(
        b1,
        confdelaunay = true, # Ensure Delaunay property
        minangle = 1, # Minimum angle (degrees)
        nosteiner = true, # Disallow new points at the boundary
        quality = false, # Don't care about grid quality, focus on Delaunay
    )
end

# ╔═╡ 38db1757-2caa-4339-9aeb-177acb9d49a5
gridplot(g_ring2)

# ╔═╡ 68f939db-6dff-4658-8995-a9f5936ef687
nondelaunay(g_ring2)

# ╔═╡ aeaac89a-5666-465b-9e38-25b980407f1a
md"""
## Step 2: Create inner part

We create the inner part using Triangulate via SimplexGridFactory. 
We insert the inner boundary of the ring as part of the geometry description.
"""

# ╔═╡ a4dba614-5955-4a02-afcd-c41160c352ea
begin
    b = SimplexGridBuilder(; Generator = Triangulate)
    bregions!(b, g_ring2, [1])
    cellregion!(b, 1)
    regionpoint!(b, 0, 0)
    g_disk = simplexgrid(b, maxvolume = hmax^2 / 2, nosteiner = true, confdelunay = true)
end

# ╔═╡ f91c6948-9228-4a8c-8aff-0c62bd5007b5
gridplot(g_disk)

# ╔═╡ dd8afc9e-f6f3-4fe3-b2da-9f257c41bc3b
nondelaunay(g_disk)

# ╔═╡ be4727bc-b1dc-4f6e-875c-3c0832d33986
md"""
## Step 3:  Glue inner part and ring
"""

# ╔═╡ 3ef68de0-1f52-46a7-8e20-a6f001060d9e
g_base = glue(g_ring2, g_disk, g1regions = [1], naive = false, strict = true)

# ╔═╡ 8922be39-5e1b-4cbe-82c9-cccd5fc34160
gridplot(g_base)

# ╔═╡ ac343bed-c301-4aa7-9a5e-331f6f01b119
nondelaunay(g_base, tol = -1.0e-17)

# ╔═╡ 436a74ae-148b-4c8f-b789-04fee27605f3
md"""
## Step 4: Tensor product with Z-coordinates
"""

# ╔═╡ 088314ef-3896-4826-b1d2-53cdd9c1f4ea
Z = range(0, len, length = nz)

# ╔═╡ f13783f1-73db-474d-82d4-da57ca91220c
g_cyl = simplexgrid(g_base, Z, top_offset = 2)

# ╔═╡ 8027688a-600f-477c-a581-6bd785a517d6
gridplot(g_cyl, Plotter = PlutoVista, zplanes = [5])

# ╔═╡ 56502247-0d07-4cba-8a77-d1e02b8195fb
nondelaunay(g_cyl, tol = 1.0e-14)

# ╔═╡ Cell order:
# ╠═784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─732e7d64-0b14-409b-9776-c222951bd79f
# ╠═84d3f058-3838-4c38-afc5-f2ea30fb9bb6
# ╟─08c71069-4ecb-450c-a905-ca6e7a62ac88
# ╠═ea4d6975-bbf9-4b72-b3df-aa4d1faea6c7
# ╠═a5120dd7-4442-4f87-90c0-bdbce949cd31
# ╠═364db500-712e-4b02-b124-1fa005e53e8e
# ╠═f4c5efff-d36a-445c-b59d-abfdc5690212
# ╟─79e6fa42-fdb2-4b99-a3ee-f5281a640fb2
# ╟─6784f6de-df9d-431b-b4b5-ca23f1bcde65
# ╠═83070579-667c-4475-9716-0980a2447fdf
# ╟─436ede6d-7c88-40f1-bf0f-b2691447314f
# ╠═2638116e-d895-4099-a6e4-2a996512cb4a
# ╠═38db1757-2caa-4339-9aeb-177acb9d49a5
# ╠═68f939db-6dff-4658-8995-a9f5936ef687
# ╟─aeaac89a-5666-465b-9e38-25b980407f1a
# ╠═a4dba614-5955-4a02-afcd-c41160c352ea
# ╠═f91c6948-9228-4a8c-8aff-0c62bd5007b5
# ╠═dd8afc9e-f6f3-4fe3-b2da-9f257c41bc3b
# ╟─be4727bc-b1dc-4f6e-875c-3c0832d33986
# ╠═3ef68de0-1f52-46a7-8e20-a6f001060d9e
# ╠═8922be39-5e1b-4cbe-82c9-cccd5fc34160
# ╠═ac343bed-c301-4aa7-9a5e-331f6f01b119
# ╟─436a74ae-148b-4c8f-b789-04fee27605f3
# ╠═088314ef-3896-4826-b1d2-53cdd9c1f4ea
# ╠═f13783f1-73db-474d-82d4-da57ca91220c
# ╠═8027688a-600f-477c-a581-6bd785a517d6
# ╠═56502247-0d07-4cba-8a77-d1e02b8195fb
