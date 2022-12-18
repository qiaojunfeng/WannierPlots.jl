using Wannier: BVectors, BVectorShells, make_supercell
using PlotlyJS
using PeriodicTable: elements

"""
    plot(bvectors::BVectors; n_k::Int=1)

Plot the b-vectors.

Plot the reciprocal lattice as arrows, the kpoints as grey dots,
and the b-vectors as spheres whose radius is proportional to the weight.
The b-vectors from the same shell have the same color.

# Arguments
- `bvectors`: the BVectors
- `n_k`: the number of repeated kpoints in each direction
"""
function plot(bvectors::BVectors; n_k::Int=1)
    traces = _plotly_lattice(bvectors.recip_lattice)

    supercell, _ = make_supercell(bvectors.kpoints, n_k)
    supercell_cart = bvectors.recip_lattice * supercell

    d = Dict(
        "type" => "scatter3d",
        "mode" => "markers",
        "x" => supercell_cart[1, :],
        "y" => supercell_cart[2, :],
        "z" => supercell_cart[3, :],
        "marker" => attr(; size=3, color="grey"),
    )
    kpts = PlotlyJS.scatter3d(d)
    push!(traces, kpts)

    mask = falses(bvectors.n_bvecs)
    ib = 1
    while ib <= bvectors.n_bvecs
        mask .= isapprox.(bvectors.weights[ib], bvectors.weights)
        vecs = bvectors.bvectors[:, mask]  # this is in Cartesian
        radius = 1e-2 * abs(bvectors.weights[ib])
        # I will just use a color from periodic table, this is simple and won't repeat
        # Hydrogen is pure white, skip it
        c = elements[ib + 1].cpk_hex
        colorscale = [[0, c], [1, c]]

        for (i, v) in enumerate(eachcol(vecs))
            idx = ib + i - 1
            hovertext = "bvec $(idx)<br>"
            hovertext *= "x: $(round(v[1]; digits=7))<br>"
            hovertext *= "y: $(round(v[2]; digits=7))<br>"
            hovertext *= "z: $(round(v[3]; digits=7))<br>"
            hovertext *= "wb: $(round(bvectors.weights[idx]; digits=7))"
            sph = _sphere(
                v,
                radius;
                colorscale=colorscale,
                hoverinfo="text",
                hovertext=hovertext,
                showscale=false,
            )
            push!(traces, sph)
        end
        ib += count(mask)
    end

    return Plot(traces, TRANSPARENT_LAYOUT)
end
