using Documenter
using Literate
using WannierPlots

makedocs(;
    sitename="WannierPlots.jl",
    authors="Junfeng Qiao and contributors.",
    modules=[WannierPlots],
    pages=[
        "Home" => "index.md",
        "API" => ["Band" => "api/band.md", "Real space" => "api/realspace.md"],
    ],
)
