function get_EvsV_curve(datafile, kinlabel::Symbol, Ecut::T) where T
    data = open(JSON3.read, datafile)
    a_list = data[kinlabel][Ecut][:LatticeConstants]
    E₀_list = data[kinlabel][Ecut][:Energies]
    a_list, E₀_list
end

function get_evx_input_file(datafile, kinlabel::Symbol, Ecut::T, filename::String) where T
    # Extract data
    @info "Reading data from file $datafile"
    a_list, E₀_list = get_EvsV_curve(datafile, kinlabel, Ecut)

    # Write input file
    file = open(filename, "w")
    for i in 1:length(a_list)
        println(file, @sprintf "%16.12f %16.12f" a_list[i] 2*E₀_list[i]) # Energies are given in Ry
    end
    close(file)
    filename
end

function read_evx_output(file)
    # Extract the raw content of the file
    xml_string = open(io->read(io, String), file)
    filename = split(file,"/")[end]
    # The ev.x program is bug and forget to close the <xml> node
    !endswith(xml_string,"</xml>\n") && (xml_string *="</xml>")

    # Extract actual data
    data = parse_xml(xml_string)
    equation_type = data["EQUATIONS_OF_STATE"]["EQUATION_TYPE"]
    Ecut = parse(Float64, split(filename, "_")[3]) # a bit dangerous
    k₀ = parse(Float64, data["EQUATIONS_PARAMETERS"]["BULK_MODULUS_KBAR"])

    # Add ?
    # E₀ = data["EQUATIONS_PARAMETERS"]["MINIMUM_ENERGY_RY"]/2 # convert from Ry to a.u.
    # V₀ = data["EQUATIONS_PARAMETERS"]["EQUILIBRIUM_VOLUME_AU_A"]

    k₀, Ecut
end

function read_evx_output_dir(dir)
    xml_files = filter(x->endswith(x, "xml"), readdir(dir))
    data = map(xml_files) do file
        read_evx_output(joinpath(dir, file))
    end
    # Sort by insceasing Ecut and return as a matrix to ease plots
    sort!(data; by=x->x[2])
    hcat([[x[i] for x in data] for i in 1:length(data[1])]...)
end
