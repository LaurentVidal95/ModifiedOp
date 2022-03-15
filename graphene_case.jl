using DelimitedFiles

include("scf_routines.jl")

# g function
c = 1/3.5; α = 1
g = x->g_cα(x, c, α)
p_g = plot_g_cα(c, α)

# scfres and band
scfres_ref = scf_graphene(;n_bands=10);
ρ_ref = scfres_ref.ρ
# # ρ_ref = regularized_density(scfres_ref.basis, ρ_ref; E_reg=20)


# K path (Γ -> K)
k_start = [0, 0, 0]; k_end = [2/3,1/3, 0]; num_k = 1000
n_bands = 5
p_path, _ = plot_info_kpath(scfres_ref.basis, k_start, k_end, num_k)
savefig(p_path, "k_path.pdf")

for Ecut in (10,20,30,40)
    @info "Ecut = $Ecut"
    # basis with regularized kinetic term
    basis = basis_PBE_graphene(RegularizedKinetic(g=g), Ecut=Ecut)
    basis_ref = basis_PBE_graphene(Kinetic(), Ecut=Ecut)

    # Adapt the reference density to the current fft_grid
    add_dim(x) = reshape(x, (size(x)...,1))
    ρ_ref_in_basis = add_dim(DFTK.interpolate_density(ρ_ref[:,:,:,1], scfres_ref.basis, basis))
    
    # Bands and derivative along path
    λ_ref, band_data_ref = bands_along_kpath(ρ_ref_in_basis, basis_ref,
                                             n_bands, k_start, k_end, num_k)
    λ, band_data = bands_along_kpath(ρ_ref_in_basis, basis, n_bands, k_start, k_end, num_k)
    ∂λ_ref = [bands_derivative(λ_ref(i), k_start, k_end, num_k) for i in 1:n_bands]
    ∂λ = [bands_derivative(λ(i), k_start, k_end, num_k) for i in 1:n_bands]
    
    p_ref = plot(1:num_k-1, ∂λ_ref[1], label="standard kinetic term")
    plot!(title="Finite difference derivative of a band w.r. to k (Graphene)"*
          "\n standard vs modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    
    p = plot(1:num_k-1, ∂λ[1], label="modified kinetic term")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    xlabel!("k path")

    p_all = plot(p_ref, p, layout=(2,1))
    plot!(size = 0.8 .*(1500,1000))

    p_bands = plot(1:num_k, λ_ref(1), label="reference band")
    plot!(1:num_k, λ(1), label="regularized band")
    plot!(size=(750,500))
    plot!(legend=:topleft)
    xlabel!("k path")
    ylabel!("ε(k)")
    
    # Save results
    output_dir = "output_2_$(Ecut)"
    !isdir(output_dir) && (mkdir(output_dir))
    writedlm(joinpath(output_dir,"energies.dat"), scfres_ref.energies.total)
    writedlm(joinpath(output_dir,"g_parameters.dat"), [c, α])    
    savefig(p_g,"$(output_dir)/g_function.pdf")
    savefig(p_all,"$(output_dir)/band_derivative_vs_kinetic_term.pdf")
    savefig(p_bands, "$(output_dir)/band_plot.pdf")
end
