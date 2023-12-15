export get_dos_plot, plot_dos

function get_dos_plot(
    x::AbstractVector,
    y::AbstractVector;
    fermi_energy=nothing,
    shift_fermi::Bool=false,
    xlabel="Energy (eV)",
    ylabel="DOS (states/eV)",
)
    traces = PlotlyJS.AbstractTrace[]

    t = PlotlyJS.scatter(;
        x,
        y,
        mode="lines",
        # line=(; color=c, kwargs...)
    )
    push!(traces, t)

    layout = PlotlyJS.Layout(;
        showlegend=false,
        xaxis=PlotlyJS.attr(;
            range=[x[1], x[end]],
            title=xlabel,
            zeroline=false,
            showgrid=true,
            showbackground=false,
            # gridcolor = KLINE_COL[],
            ticks="outside",
            showline=true,
            mirror=true,
            linecolor="black", # show axis boundary
        ),
        yaxis=PlotlyJS.attr(;
            range=[minimum(y) - 0.5, maximum(y) + 0.5],
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

    shapes = []
    # add horizontal line for y = 0
    push!(
        shapes,
        PlotlyJS.hline(
            0;
            mode="lines",
            line=PlotlyJS.attr(; dash="dash", color="grey", width=0.2),
        ),
    )
    if fermi_energy !== nothing
        if shift_fermi
            x .-= fermi_energy
            ylabel = "E - E_F (eV)"
            εF_plot = 0
        else
            εF_plot = fermi_energy
        end
        # add vertical line for Fermi to the background
        push!(
            shapes,
            PlotlyJS.vline(
                εF_plot;
                mode="lines",
                line=PlotlyJS.attr(; dash="dash", color="blue", width=0.2),
            ),
        )
    end
    PlotlyJS.relayout!(layout; shapes)

    return PlotlyJS.Plot(traces, layout)
end

"""
Plot density of states.
"""
function plot_dos(x::AbstractVector, y::AbstractVector; kwargs...)
    P = get_dos_plot(x, y; kwargs...)
    return PlotlyJS.plot(P)
end
