using PeriodicTable: elements
# prepend . due to Requires.jl
using GeometryBasics: Point, TriangleFace
using Meshing: MarchingCubes, MarchingTetrahedra, isosurface
using StatsBase: fit, Histogram
using Wannier: RGrid, origin, span_vectors

export plot_wannierfunction

# unfortunately, the PeriodicTable.jl does not provide atomic radius info,
# we can use this json from this repo
# https://github.com/AlexGustafsson/molecular-data/blob/master/json/elements.json
# TODO once there is a julia package providing atomic radius and color, switch to that one
# You need to first run
#   wget https://raw.githubusercontent.com/AlexGustafsson/molecular-data/master/json/elements.json
# elements = JSON.parsefile(joinpath(@__DIR__, "elements.json"))

"""
Guess iso from histogram.
"""
function _guess_isolevel(data)
    h = fit(Histogram, vec(data); nbins=100)
    percent = cumsum(h.weights) / length(data)
    i = findfirst(percent .>= 0.97)
    # since length(h.edges) = length(h.weights) + 1
    return h.edges[1][i + 1]
end

"""
Generate a mesh3d.
"""
function _mesh3d(
    origin::AbstractVector{T}, lattice::AbstractMatrix{T}, W::AbstractArray{T,3}, iso::T
) where {T<:Real}
    algo = MarchingCubes(; iso=iso, insidepositive=iso >= 0)
    # use marching tetrahedra with iso
    # algo = MarchingTetrahedra(; iso=iso, insidepositive=iso >= 0)
    # use Naive Surface Nets with iso
    # algo = NaiveSurfaceNets(; iso=iso, insidepositive=iso>=0)

    # generate the mesh using marching cubes
    # mc = Mesh(cube.W, algo)
    # save the file as a PLY file (change extension to save as STL, OBJ, OFF)
    # save("mc.ply", mc)

    # call isosurface to get a vector of points and vector of faces indexing to the points
    widths = [1.0, 1.0, 1.0]
    O = [0.0, 0.0, 0.0]
    vertices, faces = isosurface(W, algo; origin=O, widths=widths)

    n_vert = length(vertices)
    xyz = zeros(Float64, 3, n_vert)

    for i in 1:n_vert
        xyz[:, i] = origin + lattice * vertices[i]
    end

    n_face = length(faces)
    ijk = zeros(Int, 3, n_face)
    for i in 1:n_face
        ijk[:, i] = faces[i]
    end

    x = xyz[1, :]
    y = xyz[2, :]
    z = xyz[3, :]
    i = ijk[1, :]
    j = ijk[2, :]
    k = ijk[3, :]

    return (; x, y, z, i, j, k)
end

"""
Generate a mesh3d for plotlyjs.
"""
function _plotly_mesh3d(
    origin::AbstractVector{T},
    lattice::AbstractMatrix{T},
    W::AbstractArray{T,3},
    iso::T,
    color,
) where {T<:Real}
    x, y, z, i, j, k = _mesh3d(origin, lattice, W, iso)
    return _plotly_mesh3d(x, y, z, i, j, k, color)
end

function _plotly_mesh3d(
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    i::AbstractVector{Int},
    j::AbstractVector{Int},
    k::AbstractVector{Int},
    color,
    opacity=1.0,
) where {T<:Real}
    # plotly starts from 0
    for a in axes(i, 1)
        i[a] -= 1
        j[a] -= 1
        k[a] -= 1
    end
    t = PlotlyJS.mesh3d(;
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        color=color,
        opacity=opacity,
        lighting=attr(; ambient=0.3, diffuse=0.8, specular=0.8, roughness=0.2),
        lightposition=attr(; x=10, y=10, z=10),
    )
    return t
end

"""
Plot volumetric data with Plotly.

E.g., realspace WFs.

data: volumetric data in 3D

# Examples

1. Read xsf file
```julia
using Wannier, WannierPlots
x = read_xsf("si2_00001.xsf");
plot_wannierfunction(x.rgrid, x.W, x.primvec, inv(x.primvec) * x.atom_positions, parse.(Int, x.atoms))
```

2. Read cube file
```julia
using Wannier, WannierPlots
cube = Wannier.read_cube("Si2_valence_00001.cube");
win = Wannier.read_win("Si2_valence.win")
plot_wannierfunction(cube.rgrid, cube.W, win.unit_cell_cart, eachcol(cube.atom_positions), cube.atom_numbers)
```
"""
function plot_wannierfunction(
    rgrid::RGrid,
    W::AbstractArray{T,3},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractVector,
    atom_numbers::AbstractVector{U};
    iso::Union{T,Nothing}=nothing,
) where {T<:Real,U<:Integer}
    crystal = get_crystal_plot(lattice, atom_positions, atom_numbers)

    if isnothing(iso)
        iso = _guess_isolevel(W)
    end

    O = origin(rgrid)
    spanvec = span_vectors(rgrid)
    a, b = minimum(W), maximum(W)
    if iso >= a
        s1 = _plotly_mesh3d(O, spanvec, W, iso, "#C3423F")
    end
    if iso <= b
        s2 = _plotly_mesh3d(O, spanvec, W, -iso, "#5BC0EB")
    end

    traces = crystal.data
    push!(traces, s1)
    push!(traces, s2)

    # TODO Whatever I tried, plotlyjs refuses to set equal aspect ratio for the 3 axes :-(
    # layout = Layout(;
    #     # scene = attr(;
    #     #     aspectmode = "cube",
    #     # #     aspectratio = attr(; x=1, y=1, z=1),
    #     # ),
    #     # yaxis=attr(;
    #     #     scaleanchor = "x",
    #     #     scaleratio = 1,
    #     # ),
    #     # zaxis=attr(;
    #     #     scaleanchor = "x",
    #     #     scaleratio = 1,
    #     # ),
    # )
    return Plot(traces)#, layout)
end
