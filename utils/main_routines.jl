# Regularizing function
reg_function = y->gm(y, ha_2(0.4, -1))

"""
Compute ground state density using a big Ecut. This density is used to compute
the effective potential of H_k is every following computations.
Also generates the standard band plot k-path for the given system.
"""
function reference_data(system; k_path_res=200, n_bands=13)
    # Launch scf with standard kinetic term
    scfres_ref = system.scf()
    kpath = high_symmetry_kpath(scfres_ref.basis.model, kline_density=k_path_res)
    @info "Computing band structure"
    band_data = compute_bands(scfres_ref.basis, kpath.kcoords;
                             n_bands=n_bands, ρ = scfres_ref.ρ)
    @info "number of k-points: $(length(kpath.kcoords))"
    # Computation of k path for band plot
    (;system=system, scfres=scfres_ref, kpath=kpath, band_data=band_data)
end

"""
Compute bands and bands derivative for both standard and regularized
kinetic term given Ecut, the number of wanted bands and the regularizing function
gm.
Ecut is to be taken small in order to show irregularities.
"ref_data" is the output of the previous function "reference_data".
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

    # Compute energy bands and band derivatives along k-path
    kcoords = ref_data.kpath.kcoords
    ρ_ref = ref_data.scfres.ρ
    # # Interpolate the reference density in the current fft_grid
    # # Use only if fft_size of basis and basis_std is different from basis_ref.fft_size
    # # Less accurate than simply taking the same fft_size for every one.
    # add_dim(x) = reshape(x, (size(x)...,1))
    # ρ_ref = add_dim(DFTK.interpolate_density(ref_data.scfres.ρ[:,:,:,1],
    #                                          ref_data.scfres.basis, basis))

    # standard
    @info "Computing band structure for standard kinetic term"
    band_data_std = compute_bands(basis_std, kcoords, n_bands=n_bands, ρ=ρ_ref)
    εn_std = n->[εnk[n] for εnk in band_data_std.λ]
    ∂εn_std = [band_derivative(εn_std(n), kcoords) for n in 1:n_bands]

    # regularized term
    @info "Computing band structure for modified kinetic term"
    band_data = compute_bands(basis, kcoords, n_bands=n_bands, ρ=ρ_ref)
    εn = n->[εnk[n] for εnk in band_data.λ]
    ∂εn = [band_derivative(εn(n), kcoords) for n in 1:n_bands]

    # Return ref_data and reg_data
    (;std_data=(band_data_std, εn_std, ∂εn_std), reg_data=(band_data, εn, ∂εn))
end
