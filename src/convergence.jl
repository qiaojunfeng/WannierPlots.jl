using GLMakie
using LaTeXStrings

"""
Plot the convergence of disentanglement and maximal localization of Wannier90.

# Keyword Arguments
- `xscale`: Scale of the x-axis. Default is `identity`, use `log10` for log scale.

```julia
using WannierIO: read_wout
filename = "Si2.wout"
wout = read_wout(filename; iterations=true)
plot_wannier_convergence(wout.iterations)
```
"""
function plot_wannier_convergence(iterations::NamedTuple; xscale=identity)
    fig = Figure()
    # empty!(fig)

    # xscale = log10
    # xscale = identity
    xgridvisible = false
    ygridvisible = false

    # disentanglement
    if haskey(iterations, :disentangle)
        ax = Axis(fig[1, 1];
            title = "Disentanglement convergence",
            xlabel = "Iteration",
            ylabel = L"\Omega_{\textrm{I}}~(\AA^2)",
            xgridvisible,
            ygridvisible,
            xscale,
        )
        x = iterations.disentangle.iter
        y = iterations.disentangle.ΩI_current
        lines!(ax, x, y; label=L"\Omega_{\textrm{I}}")
        row = 2
    else
        row = 1
    end

    # max localization
    ax1 = Axis(fig[row, 1];
        title = "Maximal localization convergence",
        xlabel = "Iteration",
        ylabel = L"\Omega~(\AA^2)",
        # yticklabelcolor = :blue,
        xgridvisible,
        ygridvisible,
        xscale,
    )
    ax2 = Axis(fig[row, 1];
        ylabel = L"r~(\AA)",
        yaxisposition = :right,
        # yticklabelcolor = :red,
        xgridvisible,
        ygridvisible,
        xscale,
    )
    hidespines!(ax2)
    hidexdecorations!(ax2)

    wan = iterations.wannierize
    x = wan.iter
    if xscale == log10 && x[1] == 0
        x = x .+ 1
    end
    line_Ω = lines!(ax1, x, wan.Ωtotal; label=L"\Omega", color=:black)

    rx = [r[1] for r in wan.sum_centers]
    ry = [r[2] for r in wan.sum_centers]
    rz = [r[3] for r in wan.sum_centers]
    line_x = lines!(ax2, x, rx, label=L"r_x")
    line_y = lines!(ax2, x, ry, label=L"r_y")
    line_z = lines!(ax2, x, rz, label=L"r_z")

    lns = [line_Ω, line_x, line_y, line_z]
    labs = [l.label[] for l in lns]
    axislegend(ax1, lns, labs)

    return fig
end
