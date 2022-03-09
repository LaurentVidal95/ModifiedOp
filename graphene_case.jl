using DelimitedFiles

include("scf_routines.jl")

# g function
c = 1/3.5; α = 1
g = x->g_cα(x, c, α)
p_g = plot_g_cα(c, α)

# scfres
# scfres_ref = scf_graphene();
# scfres = scf_graphene(terms=modified_PBE_terms(g=g));

# Kpath
k_start = [0, 0, 0]; k_end = [rand()-1/2, rand()-1/2, 0]; num_k = 1000
n_bands = 5
p_path, _ = plot_info_kpath(scfres_ref.basis, k_start, k_end, num_k)

# Bands and derivative along path
# λ_ref, band_data_ref = bands_along_kpath(scfres_ref, n_bands, k_start, k_end, num_k)
# λ, band_data = bands_along_kpath(scfres, n_bands, k_start, k_end, num_k)
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

# Save results
output_dir = "output_graphene/test_1"
!isdir(output_dir) && (mkdir(output_dir))
savefig(p_path, "$(output_dir)/k_path.pdf")
savefig(p_g,"$(output_dir)/g_function.pdf")
savefig(p_all,"$(output_dir)/band_derivative_vs_kinetic_term.pdf")

writedlm(joinpath(output_dir,"energies.dat"), [scfres_ref.energies.total, scfres.energies.total])
writedlm(joinpath(output_dir,"g_parameters.dat"), [c, α])
