# Documentation

## Writing

The documentation is written in Markdown, then processed by `Documenter.jl`.

The tutorials are written in plain Julia code, then processed by `Literate.jl`,
converted to Markdown, Jupyter notebook, and Julia script,
finally the Markdown is processed by `Documenter.jl`.
Thus the tutorials are placed in the `tutorials_src` folder,
and the `Literate.jl`-generated Markdowns are placed in the `tutorials` folder.

The quickest starting point is looking at existing tutorials,
e.g. [`1-band.jl`](src/tutorials_src/1-band.jl).

## Build

Build docs locally

```shell
julia --project=. make.jl; python -m http.server --directory build
```
