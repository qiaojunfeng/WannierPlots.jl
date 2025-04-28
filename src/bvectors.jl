using Wannier: KspaceStencil, n_bvectors, make_supercell
using PeriodicTable: elements

export plot_kstencil

"""
    $(SIGNATURES)

Plot the b-vectors.

Plot the reciprocal lattice as arrows, the kpoints as grey dots,
and the b-vectors as spheres whose radius is proportional to the weight.
The b-vectors from the same shell have the same color.

# Arguments
- `kstencil`: the KspaceStencil
- `n_k`: the number of repeated kpoints in each direction
"""
function plot_kstencil(kstencil::KspaceStencil; n_k::Int=1)
    traces = _lattice(kstencil.recip_lattice; label_prefix="b")

    # This generates too much points, hard to see
    supercell, _ = make_supercell(kstencil.kpoints, n_k)
    # I just translate by half of the reciprocal lattice
    # in fractional coordinates
    # supercell = kstencil.kpoints .- Ref([0.5, 0.5, 0.5])
    # supercell, _ = make_supercell(supercell, n_k)
    supercell_cart = Ref(kstencil.recip_lattice) .* supercell

    d = Dict(
        "type" => "scatter3d",
        "mode" => "markers",
        "x" => [k[1] for k in supercell_cart],
        "y" => [k[2] for k in supercell_cart],
        "z" => [k[3] for k in supercell_cart],
        "marker" => attr(; size=3, color="grey"),
    )
    kpts = PlotlyJS.scatter3d(d)
    push!(traces, kpts)

    nbs = n_bvectors(kstencil)
    max_w = maximum(kstencil.bweights)
    mask = falses(nbs)
    ib = 1
    while ib <= nbs
        mask .= isapprox.(kstencil.bweights[ib], kstencil.bweights)
        vecs = kstencil.bvectors[mask]  # this is in Cartesian
        # radius = 1
        radius = 0.2 * (abs(kstencil.bweights[ib]) / max_w)^(1/2)
        # I will just use a color from periodic table, this is simple and won't repeat
        # Hydrogen is pure white, and the next color is also too light, skip them
        c = elements[ib + 2].cpk_hex
        colorscale = [[0, c], [1, c]]

        for (i, v) in enumerate(vecs)
            idx = ib + i - 1
            hovertext = "bvec $(idx)<br>"
            hovertext *= "x: $(round(v[1]; digits=7))<br>"
            hovertext *= "y: $(round(v[2]; digits=7))<br>"
            hovertext *= "z: $(round(v[3]; digits=7))<br>"
            hovertext *= "wb: $(round(kstencil.bweights[idx]; digits=7))"
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

    return PlotlyJS.Plot(traces, TRANSPARENT_LAYOUT)
end
