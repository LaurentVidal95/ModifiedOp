# SCF
using DFTK
using LinearAlgebra

# Plots
using Plots
using StatsPlots
using Unitful, UnitfulAtomic
using LaTeXStrings

# Misc
using StaticArrays
using JSON3
using ProgressMeter
using Printf

# 𝒢 function
using ForwardDiff

# Include utils dir
function include_dir(dir)
    dir_content = joinpath.(Ref(dir), readdir(dir))
    for thing in dir_content
        # if its a file include it
        if !isdir(thing)
            @info "included $thing"
            include(thing)
        else
            include_dir(thing)
        end
    end
    nothing    
end

include_dir("deps")
