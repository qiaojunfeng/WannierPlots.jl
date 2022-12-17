# # 3. Realspace WFs of graphene

#=
```@meta
CurrentModule = WannierPlots
```
=#

#=
In this tutorial, we will plot Fermi surface.

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`3-fermisurf.ipynb`](./3-fermisurf.ipynb)
    - Julia script: [`3-fermisurf.jl`](./3-fermisurf.jl)
=#

# ## Preparation
# Load the package
using Wannier
using WannierPlots

# Path of current tutorial
CUR_DIR = "8-fermi_surface"

bxsf = Wannier.read_bxsf("$CUR_DIR/cu.bxsf")

fig = WannierPlots.plot_fermisurf_plotly(bxsf.rgrid, bxsf.fermi_energy, bxsf.E)
fig.layout.width = 1000
fig.layout.height = 1000
fig.layout.autosize = false
fig
