Julia code to assert accuracy and regularizing properties of hamiltonian with blow up
kinetic terms.
This code is used in [1] -> [insert arXiv link]


# Requirements:
Julia 1.7 with the libraries:
- Plots, Measures, LaTeXStrings for plottings;
- DelimitedFiles for saving data;
This code is a simple top layer over a custom version of the package DFTK [insert ref] (provided in this repository).

# Usage
This code runs the computations to perform the numerical experiments from [1].
To perform the computations, first open the Julia shell with `julia --project` in your
local copy of this repository and simply call
```
include("launch_computations.jl")
```
Default system is silicon with PBE exchange correlation functional.
The code supports all other test cases defined in `utils/chemical_systems.jl`, by simply
redefining the variable `system=...` in `launch_computations.jl`.

# Contact
This is research code, not necessarily user-friendly, actively maintened or extremely robust.
If you have questions or encounter problems, contact [instert contact]

