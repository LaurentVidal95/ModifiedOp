function store_energy_curve_data(output::AbstractDict, output_file::String)
    @info "Storing data in file $(output_file)"
    # Extract computations information
    params = output["Parameters"]
    Ecut = params.Ecut
    TR = typeof(Ecut)
    blowup_rate = TR(params.blowup_rate)

    write_json(data, file) = open(io->JSON3.write(io, data, allow_inf=true), file, "w")
    parse_JSON3_key(key) = parse(TR, String(key))

    # If the file 
    if !isfile(output_file)
        data = Dict{TR, Any}()
        data[blowup_rate] = Dict{TR, Any}()
        data[blowup_rate][Ecut] = output
        write_json(data, output_file)
    else
        # the raw data has JSON3 types everywhere which makes
        # it difficult to add a new entry. To avoid types conflict
        # we copy the raw_data as a new_dictionary with better types.
        newdata = parse_nested_JSON_dict(open(JSON3.read, output_file), TR)
        ref_blowuprates = collect(keys(newdata))
        # Add new blowup rate key if needed
        if !(blowup_rate ∈ ref_blowuprates)
            newdata[blowup_rate] = Dict{TR, Any}()
        end
        newdata[blowup_rate][Ecut] = output
        write_json(newdata, output_file)
    end
    nothing
end

# Needed since its hard to add an entry to a JSON dict.
function parse_nested_JSON_dict(dict::AbstractDict, TR::DataType)
    # The sort is just to keep the natural order but doesn't matter.
    # JSON3 types first have to be converted to String before parsing...
    parse_JSON3_key(key) = parse(TR, String(key))

    new_dict = Dict{Any, Any}()
    for key in keys(dict)
        if (typeof(dict[key]) <: AbstractDict) && (key≠:Parameters)
            key_TR = parse_JSON3_key(key)
            new_dict[key_TR] = parse_nested_JSON_dict(dict[key], TR)
        else
            new_dict[key] = dict[key]
        end
    end
    new_dict
end
