"""
Compute ground state density using a big Ecut. This density is used to compute
the effective potential of H_k is every following computations.
Also generates the standard band plot k-path for the given system.
"""
function reference_data(system; k_path_res=200)
    # Launch scf with standard kinetic term
    scfres_ref = system.scf()
    # Remove extra band (added for convergence)
    n_bands = length(scfres_ref.eigenvalues[1]) - 3
    kpath = high_symmetry_kpath(scfres_ref.basis.model, kline_density=k_path_res)

    @info "Computing reference band structure"
    band_data = compute_bands(scfres_ref.basis, kpath.kcoords;
                              n_bands=n_bands, ρ = scfres_ref.ρ,
                              tol=1e-10)
    # Computation of k path for band plot
    (;system=system, scfres=scfres_ref, kpath=kpath, band_data=band_data)
end

"""
Compute bands and bands derivative for both standard and modified
kinetic term given Ecut, the number of wanted bands and the regularizing function
gm.
Ecut is to be taken small in order to show irregularities.
"ref_data" is the output of the previous function "reference_data".
"""
function bandstructure_data(Ecut::T, n_bands::Int64, blow_up_function;
                            ref_data, interval=[0.7, 0.75],
                            tol=1e-10) where {T<:Real}

    @info "Compute standard and modified kinetic term bands for low Ecut = $Ecut"
    # Compute PlaneWaveBasis for given Ecut with std and modified kinetic term.
    # The global fft_grid is the same than reference in order to retain
    # all precision on the reference density to assemble the hamiltonian blocks
    ref_fft_size = ref_data.scfres.basis.fft_size
    basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
    basis = ref_data.system.basis(ModifiedKinetic(blow_up=blow_up_function), Ecut=Ecut,
                                  fft_size=ref_fft_size)

    # Compute energy bands along reference k-path
    kcoords = ref_data.kpath.kcoords
    ρ_ref = ref_data.scfres.ρ

    band_data_std = compute_bands(basis_std, kcoords; n_bands=n_bands, ρ=ρ_ref, tol)
    band_data_mod = compute_bands(basis, kcoords; n_bands=n_bands, ρ=ρ_ref, tol)
    εn_std = n->[εnk[n] for εnk in band_data_std.λ]
    εn_mod = n->[εnk[n] for εnk in band_data_mod.λ]

    # Return ref_data and mod_data
    (;std_data=(band_data_std, εn_std), mod_data=(band_data_mod, εn_mod))
end

function extract_blow_up_rate(model)
    (model.model_name=="ModifiedKinetic") &&
        (return model.term_types[1].blow_up_function.g3.ε)
    NaN
end
extract_blow_up_rate(basis::PlaneWaveBasis) = extract_blow_up_rate(basis.model)

function focus_on_band(n, basis_in; ref_data,
                       num_k=100, path_section=ref_data.kpath.kpath[1][1:2],
                       tol=1e-4,
                       maxiter=200,
                       debug=5)
    # Compute zone to zoom on
    kpath = ref_data.kpath
    k_start_label, k_end_label = path_section
    blow_up_rate = extract_blow_up_rate(basis_in)
    @info "Focusing on band $n between $(k_start_label) and $(k_end_label) with "*
          "$(num_k) points.\n"*"Blow-up rate: $(blow_up_rate)"

    # Pre-computations
    kcoords = generate_kpath(kpath.klabels[k_start_label], kpath.klabels[k_end_label], num_k)
    ρ_ref = ref_data.scfres.ρ

    # Compute band with higher accuracy between two band diagram points
    # modified kinetic term
    band_data = compute_bands(basis_in, kcoords, n_bands=n, ρ=ρ_ref, tol=tol, maxiter=maxiter)
    εn = [εnk[n] for εnk in band_data.λ]

    # Plot finite diff derivatives
    # Put in plot directly...
    subset(tab, n) = [x for (i,x) in enumerate(tab) if rem(i,n)==0]
    tmp_εn, tmp_kcoords = subset.((εn, kcoords), Ref(debug))
    ∂εn = band_derivative(tmp_εn, tmp_kcoords)
    ∂2εn = band_derivative(∂εn, tmp_kcoords)

    (;path_section, data=(band_data, εn, ∂εn, ∂2εn))
end

# Work in progress to compare density of state
function plot_dos_perso(data)
    basis = data[1][1]
    eigenvalues = data[1][2]
    εF = DFTK.fermi_level(basis, eigenvalues)
    plot_dos(basis, eigenvalues; εF)
end

# # Interpolate the reference density in the current fft_grid
# # Use only if fft_size of basis and basis_std is different from basis_ref.fft_size
# # Less accurate than simply taking the same fft_size for every one.
# add_dim(x) = reshape(x, (size(x)...,1))
# ρ_ref = add_dim(DFTK.interpolate_density(ref_data.scfres.ρ[:,:,:,1],
#                                          ref_data.scfres.basis, basis))

function compare_dos(plot_data; ref_data)
    function compute_idos(band_data)
        basis = band_data.basis
        eigenvalues = band_data.λ        
        n_spin = basis.model.n_spin_components
        εs = range(minimum(minimum(eigenvalues)) - .5,
                   maximum(maximum(eigenvalues)) + .5, length=1000)
        Dεs = compute_dos.(εs, Ref(basis), Ref(eigenvalues))
        D = []
        for σ in 1:n_spin
            push!(D, [Dσ[σ] for Dσ in Dεs])
        end
        D
    end
    compute_idos.([plot_data.std_data[1], plot_data.mod_data[1], ref_data.band_data])
end
