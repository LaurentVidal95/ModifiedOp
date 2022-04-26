include("scf_routines.jl")

# Regularizing function
a = 0.4
ε = -1
# reg_function = y->gm(y, optimized_ha(ha_1, 1))
reg_function = y->gm(y, ha_2(a, ε))

"""
Compute ground state density using a big Ecut. This density is used to compute
the effective potential of H_k is the following computations.
"""
function reference_data(system; k_path_res=200, n_bands=13)
    # Launch scf with standard kinetic term
    scfres_ref = system.scf()
    kpath = high_symmetry_kpath(scfres_ref.basis.model, kline_density=k_path_res)
    band_data = compute_band(scfres_ref.basis, kpath.kcoords;
                             n_bands=n_bands, ρ = scfres_ref.ρ)
    @info "number of k-points: $(length(kpath.kcoords))"
    # Computation of k path for band plot
    (;system=system, scfres=scfres_ref, kpath=kpath, band_data=band_data)
end

"""
Compute bands and bands derivative for both standard and regularized
kinetic term given Ecut, the number of bands and the regularizing function
gm. "ref_data" is the output of reference_data
"""
function compute_band_and_derivatives(Ecut::T, n_bands::Int64, gm_function;
                                      ref_data) where {T<:Real}
    @info "Ecut = $Ecut"
    # Compute PlaneWaveBasis for given Ecut with std and regularized kinetic term.
    # The global fft_grid is the same than reference in order to retain
    # all precision on the reference density to assemble the hamiltonian blocks
    ref_fft_size = ref_data.scfres.basis.fft_size
    basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
    basis = ref_data.system.basis(RegularizedKinetic(g=gm_function), Ecut=Ecut,
                                  fft_size=ref_fft_size)

    # Interpolate the reference density in the current fft_grid
    # add_dim(x) = reshape(x, (size(x)...,1))
    # ρ_ref = add_dim(DFTK.interpolate_density(ref_data.scfres.ρ[:,:,:,1],
    #                                          ref_data.scfres.basis, basis))

    # Compute energy bands and band derivatives along k-path
    kcoords = ref_data.kpath.kcoords
    ρ_ref = ref_data.scfres.ρ

    # standard
    band_data_std = compute_bands(basis_std, kcoords, n_bands=n_bands, ρ=ρ_ref)
    λ_std = i->[εk[i] for εk in band_data_std.λ]
    ∂λ_std = [bands_derivative(λ_std(i), kcoords) for i in 1:n_bands]

    # regularized term
    band_data = compute_bands(basis, kcoords, n_bands=n_bands, ρ=ρ_ref)
    λ = i->[εk[i] for εk in band_data.λ]
    ∂λ = [bands_derivative(λ(i), kcoords) for i in 1:n_bands]

    # Return ref_data and reg_data
    (;ρ_ref=ρ_ref, std_data=(basis_std, band_data_std, λ_std, ∂λ_std), reg_data=(basis, band_data, λ, ∂λ))
end

function bandplots(plot_data; ref_data, savedir="")
    # Extract data
    ρ_ref = plot_data.ρ_ref
    basis_ref, band_data_std, λ_ref, ∂λ_ref = plot_data.std_data
    basis, band_data, λ, ∂λ = plot_data.reg_data

    kpath = ref_data.kpath
    kcoords = kpath.kcoords
    num_k = length(kcoords)

    # Small function that computes the fermi level and min and max eigenvalue
    # so that the shift due to the change of kinetic term does not appear.

    function fermi_and_ylims(data, kcoords)
        εF = DFTK.fermi_level(data[2].basis, data[2][2])
        ymin = findmin(data[2][2])[1][1]
        ymax = findmax(data[2][2])[1][end]
        εF, (ymin, ymax)
    end

    # Plot band diagram for standard and regularized kinetic term

    # Default plot parameters
    default(fontfamily="Computer Modern",
            linewidth=0.75, framestyle=:box,
            label=nothing, grid=:true,
            linecolor=RGB(95/255,133/255,255/255),
            margin=2mm)

    # Reference
    p_ref = DFTK.plot_band_data(ref_data.band_data;
                                ref_data.scfres.εF, kpath.klabels)
    plot!(p_ref, size=(800,500))
    ylabel!(p_ref, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f")
    title!(p_ref,"Reference band diagram")

    # Standard
    εF, (ymin, ymax) = fermi_and_ylims(plot_data.std_data, kcoords)
    p_std = DFTK.plot_band_data(band_data_std; εF, kpath.klabels)
    plot!(p_std, size=(800,500))
    ylims!(p_std, (ymin + ymin/10, ymax + ymax/10))
    xlabel!(p_std, " ")
    ylabel!(p_std, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f")
    title!(p_std,"Standard kinetic term")

    # Regularized
    εF, (ymin, ymax) = fermi_and_ylims(plot_data.reg_data, kcoords)
    p_reg = DFTK.plot_band_data(band_data; εF, kpath.klabels)
    plot!(p_reg, size=(800,500))
    ylims!(p_reg, (ymin + ymin/10, ymax + ymax/10))
    ylabel!(p_reg, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f")
    title!(p_reg,"Modified kinetic term")

    # Gather both plots
    p_bands = plot(p_std, p_reg, layout=(2,1))
    !isempty(savedir) && (savefig(p_bands, joinpath(savedir, "band_plot.pdf")))

    p_ref, p_bands
end

function plot_band_irregularity(i_band, plot_data, ref_data)
    # TODO
end
    # # Compute bands largest irregularity in order to zoom on it
    # k_irr, (k_zoom_start, k_zoom_end) = bands_irregularity(λ_ref(band_to_plot), k_coords)
    # k_coords_zoom = discretize_kpath(k_zoom_start, k_zoom_end, num_k)

    # # Compute bands around irregularity to x100 scale
    # tmp_band_data = compute_bands(basis_std, kcoords_zoom, n_bands=n_bands, ρ=ρ_ref)
    # λ_zoom_std = i->[εk[i] for εk in tmp_bands_data.λ]
    # tmp_band_data = compute_bands(basis, kcoords_zoom, n_bands=n_bands, ρ=ρ_ref)
    # λ_zoom = i->[εk[i] for εk in tmp_bands_data.λ]

    # # Plot zoom
    # p_zoom = plot(1:num_k, λ_zoom_std(band_to_plot), label="reference band")
    # plot!(1:num_k, λ_zoom(band_to_plot), label="regularized band")
    # plot!(size=(750,500))
    # plot!(legend=:topright)
    # ylabel!("ε_1(k)")
    # title!("Zoom on largest irregularity of band n°$(band_to_plot)")

    # # Plot derivatives
    # # reference
    # p_∂ref = plot(1:num_k-1, ∂λ_ref[1], label="standard kinetic term")
    # plot!(title="Finite difference derivative of the first band of Graphene along Γ → K"*
    #       "\n standard vs modified kinetic term")
    # plot!(size=(750,500))
    # plot!(legend=:topleft)
    # # regularized
    # p_∂reg = plot(1:num_k-1, ∂λ[1], label="modified kinetic term")
    # plot!(size=(750,500))
    # plot!(legend=:topleft)
    # xticks!([1,resolution,2*resolution,3*resolution], [k_start[2], k_mid[2], k_end[2]])

    # # Gather plots
    # p_∂ = plot(p_∂ref, p_∂reg, layout=(2,1))
    # plot!(size = 0.8 .*(1500,1000))

    # # Save plots and data if a dir is given
    # if !(isempty(savedir))
    #     !isdir(savedir) && (mkdir(savedir))
    #     writedlm(joinpath(savedir,"k_path.dat"), k_coords)
    #     writedlm(joinpath(savedir,"k_path_zoom.dat"), k_coords_zoom)
    #     writedlm(joinpath(savedir,"g_parameters.dat"), [c, α])
    #     savefig(p_∂,"$(savedir)/band_derivative_vs_kinetic_term.pdf")
    #     savefig(p_bands, "$(savedir)/bands_plot.pdf")
    #     savefig(p_zoom, "$(savedir)/bands_plot_zoom.pdf")
    # end

    # (p_bands, p_zoom, p_∂), (k_coords_zoom, λ_zoom_ref, λ_zoom)
#     p_bands
# end
