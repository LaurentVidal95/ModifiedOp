function FD_derivative(tab, Δ, deg=1)
    @assert (deg ≥ 0)
    (deg==0) && return tab
    ∂tab = [(tab[i+1] - tab[i-1]) / (2*Δ) for i in 2:length(tab)-1] # centered FD for O(h^2) error
    FD_derivative(∂tab, Δ, deg-1)
end

function compare_EvsV_curves(reffile, deg=0)
    data = open(JSON3.read, reffile)
    TR = Float64

    # Extract list of blowup rates and Ecuts
    blowuprates = sort(parse.( TR, String.(collect(keys(data))) ))
    parse_Ecuts(ε) = sort(parse.(TR, String.(collect(keys(data[ε])))))
    Ecut_ref = parse_Ecuts(blowuprates[1])[end]

    # Extract reference energy (standard kinetic with highest Ecut)
    _, Es_ref = get_EvsV_curve(reffile, 0.0, Ecut_ref)

    norms = Dict{TR, Any}()
    mean(tab) = sum(tab)/length(tab)
    for ε in blowuprates
        Ecuts_ε = parse_Ecuts(ε)
        norms[ε] = map(Ecuts_ε) do Ecut
            a_list, Es = get_EvsV_curve(reffile, ε, Ecut)
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

function search_plot_1(reffile, deg=0)
    res = compare_EvsV_curves(reffile_si, deg)
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

function search_plot_2(reffile, Ecut, deg=0;
                       plot_kwargs...)
    data = open(JSON3.read, reffile)
    TR = Float64

    # Extract list of blowup rates and reference Ecut
    blowuprates = sort(parse.( TR, String.(collect(keys(data))) ))
    parse_Ecuts(ε) = sort(parse.(TR, String.(collect(keys(data[ε])))))
    Ecut_ref = parse_Ecuts(blowuprates[1])[end]

    # Extract energies (standard kinetic with highest Ecut)
    as, Es_ref = get_EvsV_curve(reffile, 0.0, Ecut_ref)
    Δa = as[2] - as[1]
    
    Es = Dict{eltype(blowuprates), Any}()
    mean(tab) = sum(tab)/length(tab)
    for ε in blowuprates
        _, Es[ε] = get_EvsV_curve(reffile, ε, Ecut)
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
    dev_symbl = ["", "∂_{a}", "∂^{2}_{a}"][deg+1]
    ylabel!(p, latexstring("|$(dev_symbl)E(a) - $(dev_symbl)E_{\\mathrm{ref}}|\\quad \\mathrm{(a.u)}"))
    (deg≠0) && (title!(p, "Centered finite differences derivatives with Δa = 4e-2"))
    plot!(legendfontsize=12, legend=:topright, legendtitle="blowup rate")
    p
end


# TODO point 3)
function HF_energy_derivative(a, system, KineticTerm)
    scfres = system.scf(; a, tol=1e-10, n_bands=8)
    function HF_energy(new_a)
        basis = scfres.basis
        lattice = (new_a/a) * basis.model.lattice
        new_model = Model(basis.model; lattice, symmetries=false)
        new_basis = PlaneWaveBasis(new_model,
                                   basis.Ecut, basis.fft_size, basis.variational,
                                   basis.kcoords_global, basis.kweights_global,
                                   basis.kgrid, basis.kshift, basis.symmetries_respect_rgrid,
                                   basis.comm_kpts, basis.architecture)
        ρ = compute_density(new_basis, scfres.ψ, scfres.occupation)
        energies = energy_hamiltonian(new_basis, scfres.ψ, scfres.occupation;
                                      ρ, scfres.eigenvalues, scfres.εF).energies
        energies.total
    end
    ForwardDiff.derivative(HF_energy, a)
end

HF_energy_second_derivative(a₀, system, KineticTerm) =
    ForwardDiff.derivative(a->HF_energy_derivative(a, system, KineticTerm), a₀)
