const blowup_color_palette = [:blue, :red, :green, :brown, :orange, :black, :yellow]

function plot_energy_curves(reffile, Ecut, deg=0;
                            plot_kwargs...)
    # Extract all curves
    data = open(JSON3.read, reffile)
    TR = Float64
    blowuprates = sort(parse.( TR, String.(collect(keys(data))) ))
    
    # Plot E=f(a) for all blowup rates
    p = plot(size=1.2 .* (600, 400))
    for (i, ε) in enumerate(blowuprates)
        as, Es = read_E_vs_lattice_constant_curve(reffile, ε, Ecut)
        Δa = as[2] - as[1]
        as = as[deg+1:end-deg]
        Es = FD_derivative(Es, Δa, deg)
        plot!(as, Es; label="$ε", linecolor=blowup_color_palette[i],
              linewidth=2, plot_kwargs...)
    end
    xlabel!(p, "Lattice constant (bohr)")
    dev_symbol = ["", "∂_{a}", "∂^{2}_{a}"][deg+1]
    ylabel!(p, latexstring("$(dev_symbol)E(a)\\quad (a.u.)"))
    plot!(legendfontsize=12, legend=:topright, legendtitle="blowup rate")
    p
end

function compute_EvsV_global_errors(reffile, deg=0)
    data = open(JSON3.read, reffile)
    TR = Float64

    # Extract list of blowup rates and Ecuts
    blowuprates = sort(parse.( TR, String.(collect(keys(data))) ))
    parse_Ecuts(ε) = sort(parse.(TR, String.(collect(keys(data[ε])))))
    Ecut_ref = parse_Ecuts(blowuprates[1])[end]

    # Extract reference energy (standard kinetic with highest Ecut)
    _, Es_ref = read_E_vs_lattice_constant_curve(reffile, 0.0, Ecut_ref)

    norms = Dict{TR, Any}()
    mean(tab) = sum(tab)/length(tab)
    for ε in blowuprates
        Ecuts_ε = parse_Ecuts(ε)
        norms[ε] = map(Ecuts_ε) do Ecut
            a_list, Es = read_E_vs_lattice_constant_curve(reffile, ε, Ecut)
            Δa = a_list[2] - a_list[1]
            # Compute eventual derivatives
            Es_out = FD_derivative(Es, Δa, deg)
            Es_ref_out = FD_derivative(Es_ref, Δa, deg)
            # Compute shift if deg=0 to compensate overestimation of energies
            shift = zeros(typeof(Δa), 2)
            (deg==0) && (shift=[mean(Es), mean(Es_ref)])
            norm((Es_out .- shift[1]) .- (Es_ref_out .- shift[2]))
        end
    end

    norms
end
function plot_global_energy_landscape_error(reffile, deg=0)
    res = compute_EvsV_global_errors(reffile_si, deg)
    p = plot(size=1.2 .* (600, 400))
    for (i,ε) in enumerate([0.0, 0.5, 1.25, 1.5, 2.5])
        plot!(p, [3, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80], res[ε][1:12];
              label="$ε", linewidth=2, linecolor=[:blue, :red, :green, :brown, :orange][i],
              yscale=:log10)
    end
    xlabel!(p, "Ecut (Ha)")
    dev_symbl = ["", "∂_{a}", "∂^{2}_{a}"][deg+1]
    ylabel!(p, latexstring("||$(dev_symbl)E - $(dev_symbl)E_{\\mathrm{ref}}||\\quad \\mathrm{(a.u)}"))
    (deg≠0) && (title!(p, "Centered finite differences derivatives with Δa = 4e-2"))
    plot!(legendfontsize=12, legend=:topright, legendtitle="blowup rate")
    p
end

function plot_local_energy_landscape_error(reffile, Ecut, deg=0;
                                           plot_kwargs...)
    data = open(JSON3.read, reffile)
    TR = Float64

    # Extract list of blowup rates and reference Ecut
    blowuprates = sort(parse.( TR, String.(collect(keys(data))) ))
    parse_Ecuts(ε) = sort(parse.(TR, String.(collect(keys(data[ε])))))
    Ecut_ref = parse_Ecuts(blowuprates[1])[end]

    # Extract energies (standard kinetic with highest Ecut)
    as, Es_ref = read_E_vs_lattice_constant_curve(reffile, 0.0, Ecut_ref)
    Δa = as[2] - as[1]

    Es = Dict{eltype(blowuprates), Any}()
    mean(tab) = sum(tab)/length(tab)
    for ε in blowuprates
        _, Es[ε] = read_E_vs_lattice_constant_curve(reffile, ε, Ecut)
    end

    p = plot(size=1.2 .* (600, 400))
    for (i,ε) in enumerate(blowuprates)
        Es_ref_plot = FD_derivative(Es_ref, Δa, deg)
        Es_plot = FD_derivative(Es[ε], Δa, deg)
        plot!(as, abs.((Es_ref_plot .- mean(Es_ref_plot)) .- (Es_plot .- mean(Es_plot)));
              label="$ε", linewidth=2,
              linecolor=[:blue, :red, :green, :brown, :orange, :black, :yellow][i],
              plot_kwargs...
              )
    end
    xlabel!(p, "lattice constant (Ha)")
    dev_symbol = ["", "∂_{a}", "∂^{2}_{a}"][deg+1]
    ylabel!(p, latexstring("|$(dev_symbol)E(a) - $(dev_symbol)E_{\\mathrm{ref}}|\\quad \\mathrm{(a.u)}"))
    (deg≠0) && (title!(p, "Centered finite differences derivatives with Δa = 4e-2"))
    plot!(legendfontsize=12, legend=:topright, legendtitle="blowup rate")
    p
end
