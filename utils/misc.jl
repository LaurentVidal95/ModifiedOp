"""
Generate k-points cartesian coordinates given starting point, end point,
and the number of wanted k-points.
"""
function generate_kpath(k_start, k_end, num_k)
    [SVector{3,Float64}(kcoord) for kcoord in
               map(x-> (1-x)*k_start .+ x*k_end, LinRange(0, 1, num_k))]    
end

"""
Compute finite difference derivative of the band εn.
εn is the vector of eigenvalues for each k-points.
"""
function band_derivative(εn::Vector{T}, kcoords) where {T<:Real}
    δk = norm(kcoords[2] .- kcoords[1],1)
    [(εn[i+1]-εn[i])/δk for i in 1:length(εn)-1]
end

"""
Save the result of "focus_on_band" routine for plotting.
"""
function save_band_data(band_data, file::String)
    @assert !(isfile(file)) "$(file) already exists" 

    # Assemble dictionary
    output = Dict{String, Any}()    
    output["path_section"] = band_data.path_section
    output["εn"] = band_data.data[2]
    output["dev_εn"] = band_data.data[3]
    output["dev2_εn"] = band_data.data[4]

    # Save in file
    @info "Saving band_data in "*file*".json"
    open(io->JSON3.write(io, output, allow_inf=true), file*".json", "w")
    nothing
end

"""
Read a json file obtained with "save_band_data".
"""
function read_band_data(file)
    data_in = open(JSON3.read, file)
    (; path_section=data_in["path_section"],
     data=([], data_in["εn"], data_in["dev_εn"], data_in["dev2_εn"]))
end

DefaultInterval = [0.85, 0.90]
