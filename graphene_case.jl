using DelimitedFiles

include("scf_routines.jl")

"""
Global variables
"""
# g function
c = 1/3.5; α = 1
g = x->g_cα(x, c, α)
p_g = plot_g_cα(c, α)

# scfres and band
# scfres_ref = scf_graphene(;n_bands=10);
# ρ_ref = scfres_ref.ρ

# # Regularize density if needed
# ρ_ref = regularized_density(scfres_ref.basis, ρ_ref; E_reg=20)

# K path (Γ -> K)
k_start = ([0, 0, 0], "Γ")
k_end = ([2/3,1/3, 0], "K")
num_k = 1000

# Number of computed bands and focus on i_band
n_bands = 5; i_band = 1

function prepare_plots(Ecut)
    @info "Ecut = $Ecut"
    # basis with regularized kinetic term
    basis = basis_PBE_graphene(RegularizedKinetic(g=g), Ecut=Ecut)
    basis_ref = basis_PBE_graphene(Kinetic(), Ecut=Ecut)
    
    # Adapt the reference density to the current fft_grid
    ρ_ref_in_basis = density_in_basis(scfres_ref, basis)
    # Bands and derivative along path
    λ_ref, band_data_ref = bands_along_kpath(ρ_ref_in_basis, basis_ref,
                                             n_bands, k_start[1], k_end[1], num_k)
    λ, band_data = bands_along_kpath(ρ_ref_in_basis, basis, n_bands, k_start[1], k_end[1], num_k)
    ∂λ_ref = [bands_derivative(λ_ref(i), k_start[1], k_end[1], num_k) for i in 1:n_bands]
    ∂λ = [bands_derivative(λ(i), k_start[1], k_end[1], num_k) for i in 1:n_bands]

    # Return ref_data and reg_data
    (basis_ref, ρ_ref_in_basis, λ_ref, ∂λ_ref), (basis, λ, ∂λ)
end

function compare_bands(ref_data, reg_data; savedir="")
    basis_ref, ρ_ref_in_basis, λ_ref, ∂λ_ref = ref_data
    basis, λ, ∂λ = reg_data
    k_coords = discretize_kpath(k_start[1], k_end[1], num_k)

    # Plot bands
    p_bands = plot(1:num_k, λ_ref(i_band), label="reference band")
    plot!(1:num_k, λ(i_band), label="regularized band")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    ylabel!("ε_1(k)")
    title!("First band of graphene along k-path Γ ⟶ K\n reference term vs new kinetic term")
    xticks!([1,length(k_coords)], [k_start[2], k_end[2]])

    # Compute bands largest irregularity in order to zoom on it
    k_irr, (k_zoom_start, k_zoom_end) = bands_irregularity(λ_ref(i_band), k_coords)
    k_coords_zoom = discretize_kpath(k_zoom_start, k_zoom_end, num_k)

    # Compute bands around irregularity to x100 scale
    λ_zoom_ref, _ = bands_along_kpath(ρ_ref_in_basis, basis_ref, n_bands,
                                     k_zoom_start, k_zoom_end, num_k)
    λ_zoom, _ = bands_along_kpath(ρ_ref_in_basis, basis, n_bands, k_zoom_start, k_zoom_end, num_k)

    # Plot zoom
    p_zoom = plot(1:num_k, λ_zoom_ref(i_band), label="reference band")
    plot!(1:num_k, λ_zoom(i_band), label="regularized band")
    plot!(size=(750,500))
    plot!(legend=:topright)
    ylabel!("ε_1(k)")
    title!("Zoom x100 on irregularity")    
    xticks!([1,length(k_coords)], ["$(k_start[2]) ⟵", "⟶ $(k_end[2])"])

    # Plot derivatives
    # ref
    p_∂ref = plot(1:num_k-1, ∂λ_ref[1], label="standard kinetic term")
    plot!(title="Finite difference derivative of the first band of Graphene along Γ → K"*
          "\n standard vs modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    # reg
    p_∂reg = plot(1:num_k-1, ∂λ[1], label="modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    xticks!([1,length(k_coords)], [k_start[2], k_end[2]])
    # Gather plots
    p_∂ = plot(p_∂ref, p_∂reg, layout=(2,1))
    plot!(size = 0.8 .*(1500,1000))

    # Save if asked
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

# function compare_terms(Ecut)
#      @info "Ecut = $Ecut"
#     # basis with regularized kinetic term
#     basis = basis_PBE_graphene(RegularizedKinetic(g=g), Ecut=Ecut)
#     basis_ref = basis_PBE_graphene(Kinetic(), Ecut=Ecut)

#     # Adapt the reference density to the current fft_grid
#     ρ_ref_in_basis = density_in_basis(scfres_ref, basis)
#     # Bands and derivative along path
#     λ_ref, band_data_ref = bands_along_kpath(ρ_ref_in_basis, basis_ref,
#                                              n_bands, k_start, k_end, num_k)
#     λ, band_data = bands_along_kpath(ρ_ref_in_basis, basis, n_bands, k_start, k_end, num_k)
#     ∂λ_ref = [bands_derivative(λ_ref(i), k_start, k_end, num_k) for i in 1:n_bands]
#     ∂λ = [bands_derivative(λ(i), k_start, k_end, num_k) for i in 1:n_bands]

#     # Plot band derivative
#     p_ref = plot(1:num_k-1, ∂λ_ref[1], label="standard kinetic term")
#     plot!(title="Finite difference derivative of a band w.r. to k (Graphene)"*
#           "\n standard vs modified kinetic term")
#     plot!(size=(750,500))
#     plot!(legend=:topleft)

#     p = plot(1:num_k-1, ∂λ[1], label="modified kinetic term")
#     plot!(size=(750,500))
#     plot!(legend=:topleft)
#     xlabel!("k path")

#     p_all = plot(p_ref, p, layout=(2,1))
#     plot!(size = 0.8 .*(1500,1000))

#     # Plot bands
#     p_bands = plot(1:num_k, λ_ref(1), label="reference band")
#     plot!(1:num_k, λ(1), label="regularized band")
#     plot!(size=(750,500))
#     plot!(legend=:topleft)
#     xlabel!("k path")
#     ylabel!("ε(k)")

#     # Plot path
#     p_path, _ = plot_info_kpath(basis_ref, k_start, k_end, num_k)

#     # Save results
#     output_dir = joinpath("output/graphene", "output_$(Ecut)")
#     !isdir(output_dir) && (mkdir(output_dir))
#     writedlm(joinpath(output_dir,"energies.dat"), scfres_ref.energies.total)
#     writedlm(joinpath(output_dir,"g_parameters.dat"), [c, α])    
#     savefig(p_g,"$(output_dir)/g_function.pdf")
#     savefig(p_all,"$(output_dir)/band_derivative_vs_kinetic_term.pdf")
#     savefig(p_bands, "$(output_dir)/band_plot.pdf")
#     savefig(p_path, "$(output_dir)/k_path.pdf")
# end

#compare_terms.([10,20,30,40])
