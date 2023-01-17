using StaticArrays: SVector
using LinearAlgebra: inv
using PlotlyJS
using Brillouin: wignerseitz, KPath, reduce_to_wignerseitz

const TRANSPARENT_LAYOUT = Layout(;
    showlegend=false,
    scene=attr(;
        xaxis=attr(;
            tickvals=[],
            zeroline=false,
            showgrid=false,
            showbackground=false,
            title=attr(; text=""),
        ),
        yaxis=attr(;
            tickvals=[],
            zeroline=false,
            showgrid=false,
            showbackground=false,
            title=attr(; text=""),
        ),
        zaxis=attr(;
            tickvals=[],
            zeroline=false,
            showgrid=false,
            showbackground=false,
            title=attr(; text=""),
        ),
        aspectmode="data",
        camera=attr(; up=attr(; x=0, z=1, y=0), center=attr(; x=0, y=0, z=0)),
        dragmode="turntable",
    ),
    margin=attr(; l=0, r=0, b=0, t=0),
    autosize=true, #false,
    #width=200, height=200,
    # transparent
    plot_bgcolor="rgba(255, 255, 255, 1)",
    paper_bgcolor="rgba(255, 255, 255, 1)",
)

"""
Plot Fermi surface with Plotly.
"""
function plot_fermisurf_plotly(
    rgrid::Wannier.RGrid,
    fermi_energy::T,
    E::AbstractArray{T,4};
    kpath::Union{KPath,Nothing}=nothing,
) where {T<:Real}
    O = Wannier.origin(rgrid)

    recip_lattice = Wannier.span_vectors(rgrid)
    bz = wignerseitz([v for v in eachcol(recip_lattice)])

    if kpath === nothing
        traces = PlotlyJS.plot(bz).plot.data
    else
        # some times there are only 7 digits in bxsf, so we use a loose tolerance
        isapprox(bz.basis, kpath.basis; atol=1e-5) || error("kpath has a different reciprocal lattice")
        traces = PlotlyJS.plot(bz, kpath).plot.data
    end

    n_bands = size(E, 1)
    for i in 1:n_bands
        Ei = @view E[i, :, :, :]
        emin = minimum(Ei)
        emax = maximum(Ei)
        if emax < fermi_energy || emin > fermi_energy
            continue
        end
        x, y, z, i, j, k = _mesh3d(O, recip_lattice, Ei, fermi_energy)
        # basis is columnwise
        basis = hcat(bz.basis...)
        # move into the BZ, since the input bxsf is computed on a mesh defined by
        # three reciprocal lattice.
        x, y, z, i, j, k = _translate_fbz(basis, x, y, z, i, j, k)
        surf = _plotly_mesh3d(x, y, z, i, j, k, "")#"#C3423F")
        push!(traces, surf)
    end

    return Plot(traces, TRANSPARENT_LAYOUT)
end

"""
    _translate_fbz(O, recip_lattice, Ei, fermi_energy)

Translate eigenvalues computed on a parallelepiped mesh to the 1st BZ mesh.
"""
function _translate_fbz(
    basis::AbstractMatrix{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    i::AbstractVector{U},
    j::AbstractVector{U},
    k::AbstractVector{U},
) where {T<:Real,U<:Integer}
    # input x,y,z is in Cartesian coordinates, but I need fractional coordinates
    inv_basis = inv(basis)
    # convert to Vector of Vector since `reduce_to_wignerseitz` requires that
    basis_vec = [v for v in eachcol(basis)]

    xyz = zeros(3)  # fractional coordinates
    xyz_r = zeros(3)  # reduced back to BZ

    # check and move points to 1BZ
    # record the translations of each point
    trans = zeros(Int, 3, length(x))
    for pt in eachindex(x)
        # to fractional
        xyz .= inv_basis * [x[pt], y[pt], z[pt]]
        xyz_r .= reduce_to_wignerseitz(xyz, basis_vec)
        # to Cartesian again
        x[pt], y[pt], z[pt] = basis * xyz_r

        trans[:, pt] = round.(xyz_r - xyz)
    end

    # if the point is moved, I need to break the faces associated to this point
    # record broken faces, since most faces are not broken, the default is false
    broken = falses(length(i))
    for fc in eachindex(i)
        # three vertices of this face: i[fc], j[fc], k[fc]
        # if they are moved together, then it's fine; but if not all of them
        # are moved using the same translation, I need to break the face
        trans[:, i[fc]] == trans[:, j[fc]] == trans[:, k[fc]] || (broken[fc] = true)
    end
    # fortunately, since the bxsf stores the periodic images of
    # the left edge as the right edge, I don't need to reconnect faces

    # now remove invalidated faces
    inew = similar(i, length(i) - count(broken))
    jnew = similar(j, length(j) - count(broken))
    knew = similar(k, length(k) - count(broken))
    fcnew = 1
    for fc in eachindex(i)
        !broken[fc] || continue

        inew[fcnew] = i[fc]
        jnew[fcnew] = j[fc]
        knew[fcnew] = k[fc]
        fcnew += 1
    end

    return x, y, z, inew, jnew, knew
end
