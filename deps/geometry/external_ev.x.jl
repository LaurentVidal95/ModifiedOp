using XMLDict
using ProgressMeter

function get_evx_input_file(datafile, blowup_rate::T, Ecut::T, filename::String) where T
    # Extract data
    a_list, E₀_list = get_EvsV_curve(datafile, blowup_rate, Ecut)

    # Write input file
    file = open(filename, "w")
    for i in 1:length(a_list)
        println(file, @sprintf "%16.12f %16.12f" a_list[i] 2*E₀_list[i]) # Energies are given in Ry
    end
    close(file)
    filename
end

function run_evx_Si(datafile, blowup_rate::T, Ecut::T;
                    output_dir) where T
    # Construct evx input file
    input_file = joinpath(output_dir, "input_$(blowup_rate)_$(Ecut)_SiPBE.txt")
    output_file = joinpath(output_dir, "output_$(blowup_rate)_$(Ecut)_SiPBE.txt")
    get_evx_input_file(datafile, blowup_rate, Ecut, input_file)

    # make tmp file
    tmp_file = joinpath(output_dir, "evx_answers.txt")
    open(tmp_file, "w") do file
        println(file, "au\nfcc\n1\n$(input_file)\n$(output_file)")
    end

    # Run evx
    cmd = pipeline(`cat $(tmp_file)`, `/home/lvidal/programs/qe-7.2/bin/ev.x`)
    read(cmd, String)

    # Erase tmp files
    rm.([input_file, tmp_file])
    txt_files = joinpath.(output_dir, filter(x->endswith(x, ".txt"), readdir(output_dir)))
    rm.(txt_files)

    nothing
end

function tmp_script_si(evx_output_dir)
    for ε in [0.0, 0.5, 1.25, 1.5, 2.5]
        output_dir_ε = joinpath(evx_output_dir, "$ε")
        !(isdir(output_dir_ε)) && mkdir(output_dir_ε)
        @info "Launching ev.x for blowup rate ε"
        Ecuts = Float64.([3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80])
        for Ecut in Ecuts
            run_evx_Si(reffile_si, ε, Ecut; output_dir=output_dir_ε)
        end
    end
    nothing
end

function read_evx_xml_file(file)
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
        read_evx_xml_file(joinpath(dir, file))
    end
    # Sort by insceasing Ecut and return as a matrix to ease plots
    sort!(data; by=x->x[2])
    hcat([[x[i] for x in data] for i in 1:length(data[1])]...)
end

function plot_EvsV_curve(system, datafile::String, blowup_rates::Vector{T}, Ecut::T) where T
    data = open(JSON3.read, datafile)
    parameters = data[blowup_rates[1]][Ecut][:Parameters]

    # Extract all E = f(V) curves.
    a_list, _ = get_EvsV_curve(datafile, blowup_rates[1], Ecut)
    E₀_list = Dict{T, Any}()
    for ε in blowup_rates
        E₀_list[ε] = get_EvsV_curve(datafile, ε, Ecut)[2]
    end

    # Compute the variation of the planewave basis' size along curve
    progress = Progress(length(a_list), desc="Estimation of PW bases sizes...")
    Ma = map(a_list) do a
        basis_a = system.basis(Kinetic(); Ecut, a, parameters.kgrid)
        next!(progress)
        sum(length(G_vectors(basis_a, k)) for k in basis_a.kpoints)
    end

    # plots
    p_EvsV = plot(size=1.3 .* (600,400))
    for (i, ε) in enumerate(blowup_rates)
        plot!(p_EvsV, a_list, E₀_list[ε], label="$ε",
              linecolor=[:blue, :red, :green, :brown, :orange][i], linewidth=2)
    end

    p_Ma = plot(size=1.3 .* (600,400))
    plot!(p_Ma, a_list, Ma, linewidth=2, color=:red)
    plot(p_EvsV, p_Ma)
end

function plot_evx_outputs(dir; blowup_rates=parse.(Float64, readdir(dir)))
    # Extract all data
    data = Dict{Float64, Any}()
    for ε in blowup_rates
        data[ε] = read_evx_output_dir(joinpath(dir,"$ε"))
    end

    # Extract ref data
    k₀_ref = data[0.0][end,1] / 10 # Divide by 10 to convert to GPa
    data[0.0] = data[0.0][1:end-1, :] # ease plot by removing ref data

    # plot
    p = plot()
    for (i,ε) in enumerate(blowup_rates)
        Ecuts = data[ε][:,2]
        plot!(p, Ecuts, abs.((data[ε][:,1] ./ 10) .- k₀_ref .+ 1e-12), yscale=:log10,
        # plot!(p, Ecuts, data[ε][:,1] ./ 10,
              label="$ε",
              linecolor=[:blue, :red, :green, :brown, :orange][i], linewidth=2)
    end
    xlabel!(p, "Ecut (Ha)")
    ylabel!(p, "log₁₀(|k₀ - k₀ʳᵉᶠ|)")
    title!(p, "k₀ʳᵉᶠ ≃ $(round(k₀_ref; digits=2)) GPa")
    plot!(p, legendtitle="Blowup rates", legendtitlefontsize=14)
    plot!(xticks=[3,5,10,20,25,30,40,50,60,70,80])
    plot!(yticks=[10, 0, 0.1, 0.001, 0.0001])
    plot!(size= 1.3 .* (600,400))
    p
end
