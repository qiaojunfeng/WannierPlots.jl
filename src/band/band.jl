using PlotlyJS
using Brillouin: KPathInterpolant
using LaTeXStrings
using Wannier: get_symm_point_indices_labels, get_linear_path

export plot_band, plot_band_diff, get_band_plot, get_band_diff_plot

"""
    $(SIGNATURES)

Merge consecutive high-symmetry points.

If two high-symmetry kpoints are neighbors, merge them into one,
with label `X|Y`, where `X` and `Y` are the original labels of
the two kpoints, respectively.

# Arguments
- `symm_point_indices`: indices of high-symmetry kpoints, start from 1
- `symm_point_labels`: labels of high-symmetry kpoints
"""
function merge_consecutive_labels(
    symm_point_indices::AbstractVector{<:Integer},
    symm_point_labels::AbstractVector{<:String},
)
    idxs = copy(symm_point_indices)
    labels = copy(symm_point_labels)

    counter = 2
    for i in eachindex(symm_point_indices)[2:end]
        if symm_point_indices[i] == symm_point_indices[i - 1] + 1
            labels[counter - 1] *= "|$(symm_point_labels[i])"
        else
            idxs[counter] = symm_point_indices[i]
            labels[counter] = symm_point_labels[i]
            counter += 1
        end
    end

    return idxs[1:(counter - 1)], labels[1:(counter - 1)]
end

"""
    $(SIGNATURES)

Convert labels of high-symmetry kpoints to unicode string.

e.g., `GAMMA` -> `Γ`, `DELTA_0` -> `Δ₀`
"""
function to_unicode(labels::AbstractVector{<:AbstractString})
    label_maps = Dict(
        "GAMMA" => "Γ",
        "DELTA" => "Δ",
        "LAMBDA" => "Λ",
        "SIGMA" => "Σ",
        "0" => "₀",
        "1" => "₁",
        "2" => "₂",
        "3" => "₃",
        "4" => "₄",
        "5" => "₅",
        "6" => "₆",
        "7" => "₇",
        "8" => "₈",
        "9" => "₉",
    )
    return map(labels) do l
        if occursin("_", l)
            base, sub = split(l, "_"; limit=2)
            return get(label_maps, base, base) * get(label_maps, sub, sub)
        end
        return get(label_maps, l, l)
    end
end

function default_label_energy(shift_fermi::Bool)
    return shift_fermi ? "E - E<sub>F</sub> (eV)" : "Energy (eV)"
end

"""
    $(SIGNATURES)

Return a PlotlyJS `Plot` struct for the band structure.

# Arguments
- `x`: Vector for x axis
- `eigenvalues`: band energies, length-`n_kpoints` vector, each element is a
    length-`n_bands` vector

# Keyword Arguments
- `fermi_energy`: Fermi energy, will draw a horizontal line
- `shift_fermi`: shift the Fermi energy to 0
- `symm_point_indices`: indices of high-symmetry kpoints
- `symm_point_labels`: labels of high-symmetry kpoints
- `color`: color of the band, can be
    - a string for the color of all the eigenvalues, e.g., `"red"`
    - a Vector having the same size as `eigenvalues`, to color each individual
        eigenvalue
- `kwargs`: additional keyword arguments for PlotlyJS, the `line` attribute of `scatter`.
    For a full customization, directly manipulate the returned `Plot` struct.
    See [PlotlyJS documentation](http://juliaplots.org/PlotlyJS.jl/stable/syncplots/).
"""
function get_band_plot(
    x::AbstractVector{<:Real},
    eigenvalues::AbstractVector;
    fermi_energy::Union{Nothing,Real}=nothing,
    shift_fermi::Bool=false,
    symm_point_indices::Union{Nothing,AbstractVector}=nothing,
    symm_point_labels::Union{Nothing,AbstractVector}=nothing,
    color="black",
    ylabel=default_label_energy(shift_fermi),
    kwargs...,
)
    nkpts = length(eigenvalues)
    @assert nkpts > 0 "empty eigenvalues"
    nbands = length(eigenvalues[1])

    if symm_point_indices !== nothing
        length(symm_point_indices) == length(symm_point_labels) ||
            error("symm_idx and symm_label must have the same length")
    end
    if shift_fermi && fermi_energy === nothing
        error("shift_fermi is true, but fermi_energy is not given")
    end

    # convert to dense matrix for fast indexing, size = nbands * nkpts
    # I use captial letters for matrices
    E = reduce(hcat, eigenvalues)
    if shift_fermi
        E .-= fermi_energy
    end

    if color isa AbstractVector
        uniform_color = false
        # size = nbands * nkpts
        C = reduce(hcat, color)
        cmin, cmax = extrema(C)
        color_plot = Vector(eachrow(C))
        mode = "markers"
        color_kwargs = (; cmin, cmax, colorbar=(;), colorscale="RdBu")
    else
        uniform_color = true
        color_plot = [color for _ in 1:nbands]
        mode = "lines"
    end

    traces = PlotlyJS.AbstractTrace[]
    for (ib, e, c) in zip(1:nbands, eachrow(E), color_plot)
        text = ["ik=$(ik) iband=$(ib)" for ik in 1:nkpts]
        if uniform_color
            t = PlotlyJS.scatter(; x, y=e, mode, line=(; color=c, kwargs...), text)
        else
            t = PlotlyJS.scatter(;
                x, y=e, mode, marker=(; color=c, color_kwargs..., kwargs...), text
            )
        end
        push!(traces, t)
    end

    layout = Layout(;
        showlegend=false,
        xaxis=attr(;
            range=[x[1], x[end]],
            zeroline=false,
            showgrid=true,
            showbackground=false,
            # gridcolor = KLINE_COL[],
            ticks="outside",
            showline=true,
            mirror=true,
            linecolor="black", # show axis boundary
        ),
        yaxis=attr(;
            range=[minimum(E) - 0.5, maximum(E) + 0.5],
            title=ylabel,
            zeroline=false,
            showgrid=false,
            showbackground=false,
            ticks="outside",
            showline=true,
            mirror="all",
            linecolor="black", # show axis boundaries on all subplots
        ),
        hovermode="closest",
        autosize=true,
        # width = 480, height = 480,
        # #margin=attr(l=50, r=5, b=15, t=10),
        plot_bgcolor="rgba(0,0,0,0)",  # transparent background
        paper_bgcolor="rgba(0,0,0,0)",  # transparent background
    )

    # for storing infinite lines
    shapes = []

    if symm_point_indices !== nothing
        symm_point_labels = to_unicode(symm_point_labels)
        idxs, labels = merge_consecutive_labels(symm_point_indices, symm_point_labels)
        # add vertial lines for high-symm points to the background
        for i in idxs
            push!(shapes, vline(x[i]; mode="lines", line=attr(; color="black", width=0.2)))
        end
        # labels on x axis
        relayout!(
            layout;
            xaxis=attr(; tickmode="array", tickvals=[x[i] for i in idxs], ticktext=labels),
        )
    end

    if fermi_energy !== nothing
        if shift_fermi
            εF_plot = 0
        else
            εF_plot = fermi_energy
        end
        # add horizontal line for Fermi to the background
        push!(
            shapes,
            hline(εF_plot; mode="lines", line=attr(; dash="dash", color="deepskyblue")),
        )
    end

    relayout!(layout; shapes)
    return PlotlyJS.Plot(traces, layout)
end

function get_band_plot(kpi::KPathInterpolant, eigenvalues::AbstractVector; kwargs...)
    x = get_linear_path(kpi)
    symm_point_indices, symm_point_labels = get_symm_point_indices_labels(kpi)
    return get_band_plot(x, eigenvalues; symm_point_indices, symm_point_labels, kwargs...)
end

"""
    $(SIGNATURES)

Plot band structure.

# Arguments
- `k`: can be either
    - a Vector of Float64 for x axis value
    - a `KPathInterpolant`
- `eigenvalues`: band energies, length-`n_kpoints` vector, each element is a
    length-`n_bands` vector

# Keyword Arguments
See also the keyword arguments of [`get_band_plot`](@ref).
"""
function plot_band(k, eigenvalues::AbstractVector; kwargs...)
    P = get_band_plot(k, eigenvalues; kwargs...)
    return PlotlyJS.plot(P)
end

function get_band_diff_plot(
    k, eigenvalues_1::AbstractArray, eigenvalues_2::AbstractArray; kwargs...
)
    P1 = get_band_plot(k, eigenvalues_1; color="grey", kwargs...)
    # add legendgroup
    for (i, t) in enumerate(P1.data)
        t[:type] == "scatter" || continue
        t[:legendgroup] = "DFT"
        t[:showlegend] = (i == 1)
        t[:name] = "DFT"
        # t[:name] = "band " * string(i)
    end
    # red and slightly thinner
    P2 = get_band_plot(k, eigenvalues_2; color="red", dash="dash", width=0.9, kwargs...)
    # add legendgroup
    for (i, t) in enumerate(P2.data)
        t[:type] == "scatter" || continue
        t[:legendgroup] = "Wan"
        t[:showlegend] = (i == 1)
        t[:name] = "Wan"
        # t[:name] = "band " * string(i)
    end
    addtraces!(P1, P2.data...)

    # add legend
    P1.layout[:showlegend] = true
    P1.layout[:legend] = (;
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0.99,
        bordercolor=:lightgrey,
        borderwidth=0.8,
        bgcolor="rgba(255,255,255,0.7)",
    )
    return P1
end

"""
    $(SIGNATURES)

Compare two band structures.

`eigenvalues_1` in grey, while `eigenvalues_2` in dashed red.

# Arguments
- `k`: can be either
    - a Vector of Float64 for x axis value
    - a `KPathInterpolant`
- `eigenvalues_1`: band energies, see [`plot_band`](@ref)
- `eigenvalues_2`: band energies, see [`plot_band`](@ref)

# Keyword Arguments
See [`plot_band`](@ref)
"""
function plot_band_diff(
    k, eigenvalues_1::AbstractArray, eigenvalues_2::AbstractArray; kwargs...
)
    P = get_band_diff_plot(k, eigenvalues_1, eigenvalues_2; kwargs...)
    return PlotlyJS.plot(P)
end
