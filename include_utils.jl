using Pkg;
# Pkg.instantiate()
Pkg.activate("./");

# SCF
using DFTK
# Plots
using Plots, Measures
using Unitful, UnitfulAtomic
using DelimitedFiles
using LaTeXStrings
# Misc
using StaticArrays
using JSON3
# g function
using ForwardDiff

# Include utils dir
include.(joinpath.(Ref("utils"), readdir("utils/")))
