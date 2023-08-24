# SCF
using DFTK
using LinearAlgebra

# Plots
using Plots, Measures
using StatsPlots
using Unitful, UnitfulAtomic
using DelimitedFiles
using LaTeXStrings
# Misc
using StaticArrays
using JSON3
using ProgressMeter
using Printf

# g function
using ForwardDiff

# Include utils dir
# include("utils/test_cases.jl")
include.(joinpath.(Ref("deps"), readdir("deps/")))

