Julia code to assert accuracy and regularizing properties of hamiltonian with blow up
kinetic terms. This code is used in [1] -> [insert arXiv link]
It is a simple top layer over a custom version of the package DFTK [insert ref] (provided in this repository).

# Requirements: 
Julia 1.8.

# Installing all dependancies
Simply open a Julia shell with `julia --project` in your local copy of this repository and call
```
using Pkg; Pkg.instantiate(".")
``` 
to install all the needed dependancies.

# Usage
This code runs the computations to perform the numerical experiments from [1].

To run the computations, first open the Julia shell with `julia --project` in your
local copy of this repository and call `include("launch_computations.jl")` both to precompile
the code and to define the global parameters.

Then a simple call to the `launch_computations` function suffices. Let us do a very rough (but fast)
computation for silicon with PBE exchange correlation functional. Simply define an output directory with
`mkdir("../silicon_PBE")` and call

```
launch_computations(silicon, blowup; 
               bandplot_res=100, single_band_res=100, output_dir="../silicon_PBE", Ecut, n_bands)
```

The arguments `bandplot_res` and `single_band_res` control the number of k-points respectively in the
band diagram and in the focus on a single band. For a precise result, they should be taken above
`200` and `2000` respectively. The results of the paper are displayed for `bandplot_res=300` and 
`single_band_res=4000`. Note that the process is time consuming for these parameters and have
been run on a cluster.

The code supports graphene with PBE functional. Any other system can be added in `utils/test_cases.jl` 
by following the same syntaxe.

# Contact
This is research code, not necessarily user-friendly, actively maintened or extremely robust.
If you have questions or encounter problems, contact [instert contact]
