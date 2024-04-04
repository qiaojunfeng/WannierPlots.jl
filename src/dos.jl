export get_dos_plot, plot_dos

"""
Get a `PlotlyBase.Plot` for density of states.

# Arguments
- `x`: energies
- `y`: DOS

# Keyword Arguments
- `fermi_energy`: plot a dashed line for Fermi energy if given
- `shift_fermi`: shift Fermi energy to 0
- `xlabel`: xaxis label
- `ylabel`: yaxis label
- `cumdos`: show cumulative DOS for number of electrons
- `swap_axes`: swap x and y axes: y axis for energy and x axis for DOS
"""
function get_dos_plot(
    x::AbstractVector,
    y::AbstractVector;
    fermi_energy=nothing,
    shift_fermi::Bool=false,
    xlabel=default_label_energy(shift_fermi),
    ylabel="DOS (states/eV)",
    cumdos::Bool=false,
    swap_axes::Bool=false,
)
    traces = AbstractTrace[]

    # Fermi energy line
    if !isnothing(fermi_energy)
        if shift_fermi
            x = x .- fermi_energy
            εF = 0
        else
            εF = fermi_energy
        end
    end

    # DOS line
    color_dos = "blue"
    if cumdos
        # color of ylabel will be the same as line color
        kwargs = attr(; mode="lines", line=attr(; color=color_dos))
    else
        kwargs = attr(; mode="lines")
    end
    if swap_axes
        t = scatter(; x=y, y=x, kwargs...)
    else
        t = scatter(; x, y, kwargs...)
    end
    push!(traces, t)

    # x y axes
    kwargs = attr(;
        showlegend=false,
        hovermode="closest",
        autosize=true,
        # width = 480, height = 480,
        # #margin=attr(l=50, r=5, b=15, t=10),
        plot_bgcolor="rgba(0,0,0,0)",  # transparent background
        paper_bgcolor="rgba(0,0,0,0)",  # transparent background
    )
    xaxis = attr(;
        range=[x[1], x[end]],
        title=xlabel,
        zeroline=false,
        showgrid=true,
        showbackground=false,
        ticks="outside",
        showline=true,
        mirror=true,
        linecolor="black", # show axis boundary
    )
    if cumdos
        yaxis_title = attr(; text=ylabel, font=attr(; color=color_dos))
    else
        yaxis_title = ylabel
    end
    yaxis = attr(;
        range=[minimum(y) - 0.5, maximum(y) + 0.5],
        title=yaxis_title,
        zeroline=false,
        showgrid=false,
        showbackground=false,
        ticks="outside",
        showline=true,
        mirror="all",
        linecolor="black", # show axis boundaries on all subplots
    )
    if swap_axes
        xaxis, yaxis = yaxis, xaxis
    end
    layout = Layout(; xaxis, yaxis, kwargs...)

    # cumulative DOS
    if cumdos
        dE = x[2] - x[1]
        cdos = cumsum(y) * dE

        # axis
        color_cumdos = "red"
        kwargs = attr(;
            title=attr(; text="Number of electrons", font_color=color_cumdos),
            showline=true,
            ticks="outside",
            linecolor="black",
        )
        if swap_axes
            pop!(layout[:xaxis], :mirror, nothing)
            layout[:xaxis2] = attr(; overlaying="x", side="top", kwargs...)
        else
            pop!(layout[:yaxis], :mirror, nothing)
            layout[:yaxis2] = attr(; overlaying="y", side="right", kwargs...)
        end
        # cum DOS line
        kwargs = attr(; mode="lines", name="n_electrons", line_color=color_cumdos)
        if swap_axes
            push!(traces, scatter(; x=cdos, y=x, xaxis="x2", kwargs...))
        else
            push!(traces, scatter(; x, y=cdos, yaxis="y2", kwargs...))
        end
    end

    # horizontal/vertial lines
    shapes = []
    # add horizontal line for y = 0
    kwargs = attr(;
        mode="lines", line=PlotlyJS.attr(; dash="dash", color="grey", width=0.2)
    )
    if swap_axes
        push!(shapes, vline(0; kwargs...))
    else
        push!(shapes, hline(0; kwargs...))
    end
    # add vertical line for Fermi to the background
    if !isnothing(fermi_energy)
        kwargs = attr(;
            name="Fermi energy",
            mode="lines",
            line=attr(; dash="dash", color="black", width=0.5),
        )
        if swap_axes
            push!(shapes, hline(εF; kwargs...))
        else
            push!(shapes, vline(εF; kwargs...))
        end
    end
    relayout!(layout; shapes)

    return Plot(traces, layout)
end

"""
Plot density of states.
"""
function plot_dos(x::AbstractVector, y::AbstractVector; kwargs...)
    P = get_dos_plot(x, y; kwargs...)
    return plot(P)
end
