Julia code to assert accuracy and regularizing properties of modified hamiltonian with blow up
kinetic terms, as studied in [CHV][^1]. It is a simple top layer over the plane-wave density functional theory
package DFTK.jl[^2].

# Requirements: 
Julia 1.8. and above.

# Installing all dependancies
Open a Julia shell with `julia --project` in your local copy of this repository and call
```
using Pkg; Pkg.instantiate(".")
``` 
to install all the needed dependancies.

# Usage
This code runs the computations to perform the numerical experiments from [CHV][^1] section 6.

To run the computations, first open the Julia shell with `julia --project` in your
local copy of this repository and call `include("launch_computations.jl")` both to precompile
the code and to define the global parameters. The default parameters corresponds to the
parameters used in [CHV][^1].

Then a simple call to the `launch_computations` function suffices to perform all numerical tests.
Let us do a very rough (but fast) computation for face centered cubic crystaline silicon with PBE exchange correlation functional. 
Simply define an output directory with
`mkdir("../silicon_PBE")` and call

```
launch_computations(silicon, blowup; 
               bandplot_res=100, single_band_res=100, output_dir="../silicon_PBE", Ecut, n_bands)
```

The arguments `bandplot_res` and `single_band_res` control the number of k-points respectively in the
full band diagram and in the focus on a single band. The code uses basic finite-difference derivatives,
which imposesto fix the parameters above `200` and `2000` respectively, to avoid numerical artifacts.
The results of the paper are displayed for `bandplot_res=300` and `single_band_res=4000`.
Note that the process is time consuming for these parameters and have been run on a cluster.

The code also supports graphene with PBE functional. Any other system can be added in `utils/test_cases.jl` 
by following the same syntaxe.

# Contact
This is research code, not necessarily user-friendly, actively maintened or extremely robust.
If you have questions or encounter problems, contact us at: Laurent(dot)vidal(at)enpc(dot)fr

[^1]: E. Cancès, M. Hassan, L. Vidal: Modified-operator method for the calculation of band diagrams of crystalline materials (arXiv:2210.00442)
[^2]: https://github.com/JuliaMolSim/DFTK.jl
