using DelimitedFiles

include("scf_routines.jl")

"""
Global variables
"""
function reference_computations(system;
                                n_bands=10,
                                k_path_res=200)
    # Launch scf with standard kinetic term
    scfres_ref = system.scf(; n_bands=n_bands)

    # Computation of k path for band plot
    recip_lat = scfres_ref.basis.model.recip_lattice
    
end

# g function
ha = y->ha_2(y, 0.4, -2)
# ha = optimized_ha(ha_1, 0.05)
g = x->gm(x, ha)
p_g = plot_gm(ha)

# scfres and band
scfres_ref = scf_graphene(;n_bands=10);
recip_lat = scfres_ref.basis.model.recip_lattice
ρ_ref = scfres_ref.ρ

# # Regularize density if needed
# ρ_ref = regularized_density(scfres_ref.basis, ρ_ref; E_reg=20)

# Band diagram path provided by DFTK
resolution=200
kpath = high_symmetry_kpath(scfres_ref.model, kline_density=resolution)
kcoords = kpath.kcoords

# Number of computed bands and single band to focus on
n_bands = 8; band_2_plot = 1

"""
Compute reference data and data for the regularized kinetic term given Ecut.
"""
function prepare_plots(Ecut)
    @info "Ecut = $Ecut"
    # basis with regularized kinetic term
    basis = basis_PBE_graphene(RegularizedKinetic(g=g), Ecut=Ecut)
    basis_ref = basis_PBE_graphene(Kinetic(), Ecut=Ecut)

    # Adapt the reference density to the current fft_grid
    ρ_ref_in_basis = density_in_basis(scfres_ref, basis)
    # Bands and derivative along path
    λ_ref, band_data_ref = bands_along_kpath(ρ_ref_in_basis, basis_ref, n_bands, k_coords)
    λ, band_data = bands_along_kpath(ρ_ref_in_basis, basis, n_bands, k_coords)
    ∂λ_ref = [bands_derivative(λ_ref(i), k_coords) for i in 1:n_bands]
    ∂λ = [bands_derivative(λ(i), k_coords) for i in 1:n_bands]

    # Return ref_data and reg_data
    (;ref_data=(basis_ref, ρ_ref_in_basis, λ_ref, ∂λ_ref), reg_data=(basis, λ, ∂λ))
end

function compare_bands(ploting_data; savedir="")
    basis_ref, ρ_ref_in_basis, λ_ref, ∂λ_ref = ploting_data.ref_data
    basis, λ, ∂λ = ploting_data.reg_data
    num_k = length(k_coords)

    # Plot bands
    p_bands = plot(1:num_k, λ_ref(band_2_plot), label="reference band")
    plot!(1:num_k, λ(band_2_plot), label="regularized band")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    ylabel!("ε_1(k)")
    title!("First band of graphene along k-path Γ ⟶ M ⟶ K ⟶ Γ\n"*
           "reference term vs new kinetic term")
    xticks!([1,resolution,2*resolution,3*resolution], [k_start[2], k_mid[2], k_end[2]])

    # Compute bands largest irregularity in order to zoom on it
    k_irr, (k_zoom_start, k_zoom_end) = bands_irregularity(λ_ref(band_2_plot), k_coords)
    k_coords_zoom = discretize_kpath(k_zoom_start, k_zoom_end, num_k)

    # Compute bands around irregularity to x100 scale
    λ_zoom_ref, _ = bands_along_kpath(ρ_ref_in_basis, basis_ref, n_bands,
                                     k_zoom_start, k_zoom_end, num_k)
    λ_zoom, _ = bands_along_kpath(ρ_ref_in_basis, basis, n_bands,
                                  k_zoom_start, k_zoom_end, num_k)

    # Plot zoom
    p_zoom = plot(1:num_k, λ_zoom_ref(band_2_plot), label="reference band")
    plot!(1:num_k, λ_zoom(band_2_plot), label="regularized band")
    plot!(size=(750,500))
    plot!(legend=:topright)
    ylabel!("ε_1(k)")
    title!("Zoom x100 on irregularity")
    xticks!([1,num_k], ["$(k_start[2]) ⟵", "⟶ $(k_start[2])"])

    # Plot derivatives
    # reference
    p_∂ref = plot(1:num_k-1, ∂λ_ref[1], label="standard kinetic term")
    plot!(title="Finite difference derivative of the first band of Graphene along Γ → K"*
          "\n standard vs modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    # regularized
    p_∂reg = plot(1:num_k-1, ∂λ[1], label="modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    xticks!([1,resolution,2*resolution,3*resolution], [k_start[2], k_mid[2], k_end[2]])

    # Gather plots
    p_∂ = plot(p_∂ref, p_∂reg, layout=(2,1))
    plot!(size = 0.8 .*(1500,1000))

    # Save plots and data if a dir is given
    if !(isempty(savedir))
        !isdir(savedir) && (mkdir(savedir))
        writedlm(joinpath(savedir,"k_path.dat"), k_coords)
        writedlm(joinpath(savedir,"k_path_zoom.dat"), k_coords_zoom)
        writedlm(joinpath(savedir,"g_parameters.dat"), [c, α])
        savefig(p_∂,"$(savedir)/band_derivative_vs_kinetic_term.pdf")
        savefig(p_bands, "$(savedir)/bands_plot.pdf")
        savefig(p_zoom, "$(savedir)/bands_plot_zoom.pdf")
    end

    (p_bands, p_zoom, p_∂), (k_coords_zoom, λ_zoom_ref, λ_zoom)
end
