function bands_along_kpath(scfres, n_bands, kcoords)
    # Compute bands on given kcoords
    basis = basis_given_kcoords(scfres.basis, kcoords)
    bands_data = compute_bands(basis, kcoords, n_bands=n_bands, ρ=scfres.ρ)
    i -> [Ek[i] for Ek in bands_data.λ], bands_data
end
bands_along_kpath(scfres, n_bands, k_start, k_end, num_k) =
    compute_bands_graphene(scfres, n_bands, discretize_k_path(k_start, k_end, num_k))

function bands_derivative(εn, kcoords)
    @show δk = norm(kcoords[2] .- kcoords[1])
    [(εn[i+1]-εn[i])/δk for i in 1:length(εn)-1]
end
