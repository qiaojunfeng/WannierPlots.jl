using LinearAlgebra
using PeriodicTable: elements
# prepend . due to Requires.jl
using GeometryBasics: Point, TriangleFace
using Meshing: MarchingCubes, MarchingTetrahedra, isosurface
using StatsBase: fit, Histogram
using PlotlyJS

# unfortunately, the PeriodicTable.jl does not provide atomic radius info,
# we can use this json from this repo
# https://github.com/AlexGustafsson/molecular-data/blob/master/json/elements.json
# TODO once there is a julia package providing atomic radius and color, switch to that one
# You need to first run
#   wget https://raw.githubusercontent.com/AlexGustafsson/molecular-data/master/json/elements.json
# elements = JSON.parsefile(joinpath(@__DIR__, "elements.json"))

"""
A plotly sphere.
"""
function _sphere(r₀, r; kwargs...)
    N = 32
    u = range(0, 2π, N)
    v = range(0, π, N)
    x = r₀[1] .+ r * cos.(u) * sin.(v)'
    y = r₀[2] .+ r * sin.(u) * sin.(v)'
    z = r₀[3] .+ r * repeat(cos.(v)'; outer=[N, 1])

    return PlotlyJS.surface(; x=x, y=y, z=z, kwargs...)
end

function _lattice(lattice::AbstractMatrix)
    # trace to loop through all the lattice parallelepipede
    xyz =
        [
            0 0 0
            1 0 0
            1 1 0
            0 1 0
            0 0 0
            0 0 1
            1 0 1
            1 0 0
            1 1 0
            1 1 1
            1 0 1
            1 1 1
            0 1 1
            0 1 0
            0 1 1
            0 0 1
        ]'
    return lattice * xyz
end

function _atoms(
    lattice::AbstractMatrix, atom_positions::AbstractMatrix, atom_numbers::AbstractVector
)
    xyz = []
    radius = []
    color = []
    label = []
    for i in axes(atom_positions, 2)
        # TODO these are from json file
        # ele = elements[atom_numbers[i]]
        # s = ele["symbol"]
        # c = ele["cpkHexColor"]
        # r = ele["radius"] / elements[1]["radius"] / 10  # normalized by Hydrogen radius
        # these from PeriodicTable
        ele = elements[atom_numbers[i]]
        s = ele.symbol
        push!(label, s)
        c = ele.cpk_hex
        push!(color, c)
        # https://github.com/JuliaPhysics/PeriodicTable.jl/issues/34
        r = 0.5  # currently no radius in PeriodicTable
        push!(radius, r)
        pos = lattice * atom_positions[:, i]
        push!(xyz, pos)
    end
    return xyz, radius, color, label
end

"""
Plotly lattice

lattice: each column is a lattice vector
origin: the overall shift of the structure, in cartesian
"""
function _plotly_lattice(lattice::AbstractMatrix, origin::AbstractVector=[0, 0, 0])
    # lattice
    xyz = _lattice(lattice)
    x = xyz[1, :] .+ origin[1]
    y = xyz[2, :] .+ origin[2]
    z = xyz[3, :] .+ origin[3]

    d = Dict(
        "mode" => "lines",
        "x" => x,
        "y" => y,
        "z" => z,
        "line" => Dict("width" => 6, "color" => "black"),
    )
    lat = PlotlyJS.scatter3d(d)

    arrow = Dict(
        "x" => lattice[1, :] .+ origin[1],
        "y" => lattice[2, :] .+ origin[2],
        "z" => lattice[3, :] .+ origin[3],
        "u" => 1 .* lattice[1, :],
        "v" => 1 .* lattice[2, :],
        "w" => 1 .* lattice[3, :],
        "anchor" => "tip", # make cone tip be at endpoint
        "sizemode" => "absolute",
        "sizeref" => 0.5,
        "hoverinfo" => ["text"],
        "text" => ["a1", "a2", "a3"],
        "colorscale" => [[0, "#FDE74C"], [1, "#FDE74C"]], # color all cones yellow
        "showscale" => false,
    )
    axs = PlotlyJS.cone(arrow)

    return [lat, axs]
end

"""
Plotly lattice and atoms

lattice: each column is a lattice vector
atom_positions: each column is an atomic position, in fractional coordinates
atom_numbers: atomic number of each atom
origin: the overall shift of the structure, in cartesian
"""
function _plotly_structure(
    lattice::AbstractMatrix,
    atom_positions::AbstractMatrix,
    atom_numbers::AbstractVector;
    origin::AbstractVector=[0, 0, 0],
)
    traces = _plotly_lattice(lattice, origin)

    # atoms
    # d = Dict(
    #     "mode" => "markers",
    #     "x" => pos[1, :],
    #     "y" => pos[2, :],
    #     "z" => pos[3, :],
    #     "marker" => Dict(
    #         "size" => 5,
    #         "color" => "orange",
    #         # colorscale => "Greens",
    #         # cmin => -20,
    #         # cmax => 50
    #     )
    # )
    # atoms = scatter3d(d)
    atoms = []
    xyz, radius, color, label = _atoms(lattice, atom_positions, atom_numbers)
    for i in axes(xyz, 1)
        s = label[i]
        c = color[i]
        r = radius[i]
        colorscale = [[0, c], [1, c]]
        pos = xyz[i]
        sph = _sphere(pos, r; colorscale=colorscale, text=s, showscale=false)
        push!(atoms, sph)
    end

    append!(traces, atoms)

    return traces
end

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
"""
function plot_wf_plotly(
    rgrid::Wannier.RGrid,
    W::AbstractArray{T,3},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{U};
    iso::Union{T,Nothing}=nothing,
) where {T<:Real,U<:Integer}
    structure = _plotly_structure(lattice, atom_positions, atom_numbers)

    if isnothing(iso)
        iso = _guess_isolevel(W)
    end

    O = Wannier.origin(rgrid)
    spanvec = Wannier.span_vectors(rgrid)
    a, b = minimum(W), maximum(W)
    if iso >= a
        s1 = _plotly_mesh3d(O, spanvec, W, iso, "#C3423F")
    end
    if iso <= b
        s2 = _plotly_mesh3d(O, spanvec, W, -iso, "#5BC0EB")
    end

    traces = [structure..., s1, s2]

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

# x = read_xsf("si2_00001.xsf");
# WannierPlots.plot_wf_plotly(x.rgrid, x.W, x.primvec, inv(x.primvec) * x.atom_positions, parse.(Int, x.atoms))
