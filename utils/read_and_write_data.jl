using JSON3

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

function read_band_data(file)
    data_in = open(JSON3.read, file)
    (; path_section=data_in["path_section"],
     data=([], data_in["εn"], data_in["dev_εn"], data_in["dev2_εn"]))
end
