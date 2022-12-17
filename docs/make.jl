using Documenter
using Literate
using WannierPlots

# I put tutorials in `tutorials_src`, and let Literate.jl
# generates markdown from them. Then Documenter.jl processes
# the markdown and renders HTML. In this way, the url of the
# tutorial pages are `tutorials/<tutorial_name>` instead of
# `tutorials_src/<tutorial_name>`.
TUTORIALS_SRCDIR = joinpath(@__DIR__, "src/tutorials_src")
TUTORIALS_OUTDIR = joinpath(@__DIR__, "src/tutorials")
TUTORIALS_BUILDDIR = joinpath(@__DIR__, "build/tutorials")

# Copied from
# https://github.com/thchr/Brillouin.jl/blob/fad88c5b6965fe4bd59e725ea60655348d36ce0f/docs/make.jl#L4
# ---------------------------------------------------------------------------------------- #
# make PlotlyJS plots showable in ```@example ``` blocks, following the approach suggested
# in https://github.com/fredrikekre/Literate.jl/issues/126
using PlotlyJS
struct HTMLPlot
    p
    h::Int # desired display height in pixels
end
HTMLPlot(p) = HTMLPlot(p, 400)
const ROOT_DIR = joinpath(@__DIR__, "build/tutorials")
const PLOT_DIR = joinpath(ROOT_DIR, "plots")
function Base.show(io::IO, ::MIME"text/html", p::HTMLPlot)
    mkpath(PLOT_DIR)
    path = joinpath(PLOT_DIR, string(hash(p) % UInt32, ".html"))
    PlotlyJS.savefig(p.p, path; format="html")
    return print(
        io,
        "<object type=\"text/html\" data=\"../$(relpath(path, ROOT_DIR))\" style=\"width:100%;height:$(p.h)px;\"></object>",
    )
end
# ---------------------------------------------------------------------------------------- #

for md in readdir(TUTORIALS_SRCDIR)
    file = joinpath(TUTORIALS_SRCDIR, md)

    Literate.markdown(file, TUTORIALS_OUTDIR)

    # the notebook needs to be executed in the correct path to read `amn` etc. files,
    # however, Literate.jl will execute the notebook in the output dir,
    # so I need to first output in workdir, then move to build dir
    workdir = joinpath(@__DIR__, "../tutorials/tutorials/")
    Literate.notebook(file, workdir)
    f = splitext(md)[1]
    mv(joinpath(workdir, "$f.ipynb"), joinpath(TUTORIALS_OUTDIR, "$f.ipynb"); force=true)

    Literate.script(file, TUTORIALS_OUTDIR)
end

makedocs(;
    sitename="WannierPlots.jl",
    authors="Junfeng Qiao and contributors.",
    modules=[WannierPlots],
    pages=[
        "Home" => "index.md",
        # the tutorials will be processed by Literate
        "Tutorial" => [
            "Band structure" => "tutorials/1-band.md",
            "Real space WFs" => "tutorials/2-wf.md",
            # "Fermi surface" => "tutorials/3-fermisurf.md",
        ],
        "API" => [
            "Band" => "api/band.md",
            "Real space" => "api/realspace.md",
            # "Fermi surface" => "api/fermisurf.md",
        ],
    ],
)
