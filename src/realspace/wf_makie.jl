using LinearAlgebra
using GeometryBasics
import Makie as MK
# if you need a standalone window, call `using GLMakie` before `using WannierPlots`,
# if you need a browser tab, call `using WGLMakie` before `using WannierPlots`.

export plot_wf

"""
    plot_wf(rgrid, W::Array{Float,3}, lattice::Matrix{Float},
        atom_positions::Matrix{Float}, atom_numbers::Vector{Int};
        iso::=nothing)

Plot volumetric data with Makie, e.g., real space WF.

# Arguments
- `rgrid`: `Wannier.RGrid`
- `W`: volumetric data in 3D
- `lattice`: columns are lattice vectors
- `atom_positions`: columns are atom positions in fractional coordinates
- `atom_numbers`: atomic numbers

# Keyword Arguments
- `iso`: isosurface value, if not given, auto calculate.
"""
function plot_wf(
    rgrid::Wannier.RGrid,
    W::AbstractArray{T,3},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{U};
    iso::Union{T,Nothing}=nothing,
) where {T<:Real,U<:Integer}
    # Somehow if I create a scene, then in browser the plot disappear after
    # mouse rotation, after refresh it appears again.
    # In the end I need to return the figure instead of scene
    # scene = Makie.Scene(resolution = (1080, 900))
    eyeposition = lattice[:, 1] * 3 + lattice[:, 2] * 1 + lattice[:, 3] * 0.5
    # eyeposition .*= 100
    # println(eyeposition)
    kwargs = (;
        projectiontype="Orthographic",
        eyeposition=eyeposition,
        lookat=Vec3f(0, 0, 0),
        upvector=Vec3f(0, 0, 1),
        # fov = 180.0,
        near=0,
        far=50000,  # do not clip
    )
    # cam3d_cad!(scene; kwargs...)

    # lattice
    xyz = _lattice(lattice)
    # lines!(scene, xyz[1, :], xyz[2, :], xyz[3, :])
    MK.lines(xyz[1, :], xyz[2, :], xyz[3, :])

    fig = MK.current_figure()
    # scene = fig.scene
    # println(scene)
    MK.cam3d_cad!(fig.scene; kwargs...)

    # lattice vectors
    # https://github.com/JuliaPlots/Makie.jl/issues/1206
    head_len = 0.4
    ps = [MK.Point3f(0, 0, 0) for i in 1:3]
    ns = [MK.Point3f(a * (norm(a) - head_len) / norm(a)) for a in eachcol(lattice)]
    color = [:red, :green, :blue]
    # arrows!(scene,
    MK.arrows!(
        ps,
        ns;
        fxaa=true, # turn on anti-aliasing
        linecolor=color,
        arrowcolor=color,
        # linewidth = 0.1,
        arrowsize=Vec3f(0.3, 0.3, head_len),
        normalize=false,
    )

    # atoms
    xyz, radius, color, label = _atoms(lattice, atom_positions, atom_numbers)
    xs = [a[1] for a in xyz]
    ys = [a[2] for a in xyz]
    zs = [a[3] for a in xyz]
    # meshscatter!(scene, xs, ys, zs, markersize = radius, color = color)
    MK.meshscatter!(xs, ys, zs; markersize=radius, color=color)

    if isnothing(iso)
        iso = _guess_isolevel(W)
    end
    a, b = minimum(W), maximum(W)

    O = Wannier.origin(rgrid)
    V = Wannier.span_vectors(rgrid)
    kwargs = (;
        shading=true,
        fxaa=true, # turn on anti-aliasing
        # base light of the plot only illuminates red colors
        ambient=Vec3f(0.3, 0.3, 0.3),
        # light from source (sphere) illuminates yellow colors
        diffuse=Vec3f(0.4, 0.4, 0.4),
        # reflections illuminate blue colors
        specular=Vec3f(1.0, 1.0, 1.0),
        # Reflections are sharp
        shininess=128.0f0,
        # lightposition = Vec3f(0.5, 0.5, 3),
    )
    if iso >= a
        x, y, z, i, j, k = _mesh3d(O, V, W, iso)
        vertices = hcat(x, y, z)
        faces = hcat(i, j, k)
        # Makie.mesh!(scene, vertices, faces; color="#FE4A49", kwargs...)
        MK.mesh!(vertices, faces; color="#FE4A49", kwargs...)
    end
    if iso <= b
        x, y, z, i, j, k = _mesh3d(O, V, W, -iso)
        vertices = hcat(x, y, z)
        faces = hcat(i, j, k)
        # Makie.mesh!(scene, vertices, faces; color="#1BE7FF", kwargs...)
        MK.mesh!(vertices, faces; color="#1BE7FF", kwargs...)
    end

    # return scene
    return fig
end
