# Tutorial on Fitting Populations to Catalogs with Selection Effects

This was a notebook for a tutorial that I gave on fitting populations to catalogs with selection effects for our GW group meeting on 2025-05-14.  The tutorial worked up to the result in the [Callister, et al. (2021)]() "Who Ordered That?" paper that found a correlation between the mass ratio and the mean effective spin ("chi effective") in the LIGO catalog.  The tutorial is in the notebook `notebooks/SelectionTutorial.ipynb` and is written in Julia.  TODO: link to run in Colab?

The canonical way to install [Julia](http://julialang.org) is via [juliaup](https://julialang.org/install/).  Once you have it installed, and have cloned this repo, fire up a julia interpreter (via the command line, in VSCode, or whatever), and install the dependencies for this project via:

```julia
> using Pkg # Activate the package model that manages dependencies
> Pkg.activate("path/to/this/repo") # Activate the package for this repo
> Pkg.instantiate() # Install the dependencies for this repo (you only have to run this the first time)
```

Now the dependencies are installed, and you can run the notebook.  If you are using VSCode, you can just open the notebook, and it will automatically use a Julia kernel; if running from the command line, you will have to install IJulia run the notebook manually:

```julia
> Pkg.add("IJulia") # Install the IJulia package
> using IJulia # Load the IJulia package
> IJulia.jupyterlab() # Start the JupyterLab server
```