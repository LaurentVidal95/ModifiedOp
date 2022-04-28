using Pkg;
Pkg.activate("./");

using DFTK
using Plots, Measures
using Unitful, UnitfulAtomic
using DelimitedFiles
using LaTeXStrings

# Include utils dir
include.(joinpath.(Ref("utils"), readdir("utils/")))
