module WannierPlots

using Reexport

@reexport using PlotlyJS
# explicitly export `plot` since both WGLMakie and PlotlyJS export "plot"
@reexport using PlotlyJS: plot
# @reexport using Makie

include("band/band.jl")
include("band/diff.jl")

include("realspace/wf_plotlyjs.jl")
include("realspace/wf_makie.jl")

end
