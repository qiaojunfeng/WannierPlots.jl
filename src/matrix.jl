using Makie

"""
Display a matrix

```julia
M = rand(3, 4)
showmat(M)
showmat(symlog(.1).(M))
```
"""
function showmat(m::AbstractMatrix{T}) where {T<:Real}
    mt = transpose(m)
    nr, nc = size(m)
    xs = 1:nc
    ys = 1:nr
    fig = Figure()
    ax = Axis(fig[1,1])
    ax.xticks = (xs, string.(xs))
    ax.yticks = (ys, string.(ys))
    # Requires python, heavy dependencies
    # using PerceptualColourMaps
    # https://colorcet.com/gallery.html#linear
    # colormap = cmap("D01")  # 0 -> white
    colormap = :vik
    mc = maximum(abs.(m))
    colorrange = (-mc, mc)
    hm = heatmap!(ax, mt; colormap, colorrange)
    ax.yreversed = true
    Colorbar(fig[1,2], hm)
    fig
end

function symlog(a)
    f(x) = if x > a
        log(x/a) + 1
    elseif x < -a
        -(log(-x/a) + 1)
    else
        x/a
    end
end
