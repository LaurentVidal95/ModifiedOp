"""
Compute ground state density using a big Ecut. This density is used to compute
the effective potential of H_k is every following computations.
Also generates the standard band plot k-path for the given system.
"""
function reference_data(system; k_path_res=200, n_bands=8)
    # Launch scf with standard kinetic term
    scfres_ref = system.scf(; n_bands)

    # Remove extra band (added for convergence)
    n_bands = scfres_ref.n_bands_converge
    kpath = irrfbz_path(scfres_ref.basis.model)
    kcoords = DFTK.Brillouin.interpolate(kpath, density=k_path_res)

    @info "Computing reference band structure"
    band_data = compute_bands(scfres_ref.basis, kcoords;
                              n_bands, ρ = scfres_ref.ρ, 
                              tol=1e-10)
    # Computation of k path for band plot
    (;system, scfres=scfres_ref, kpath, kcoords, band_data)
end

"""
Compute bands and bands derivative for both standard and modified
kinetic term given Ecut, the number of wanted bands and the regularizing function
gm.
Ecut is to be taken small in order to show irregularities.
"ref_data" is the output of the previous function "reference_data".
"""
function bandstructure_data(Ecut::T, n_bands::Int64, blowup;
                            ref_data, interval=DefaultInterval,
                            tol=1e-10) where {T<:Real}

    @info "Compute standard and modified kinetic term bands for low Ecut = $Ecut"
    # Compute PlaneWaveBasis for given Ecut with std and modified kinetic term.
    # The global fft_grid is the same than reference in order to retain
    # all precision on the reference density to assemble the hamiltonian blocks
    ref_fft_size = ref_data.scfres.basis.fft_size
    basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
    basis = ref_data.system.basis(Kinetic(;blowup), Ecut=Ecut, fft_size=ref_fft_size)

    # Compute energy bands along reference k-path
    kcoords = ref_data.kcoords
    ρ_ref = ref_data.scfres.ρ

    band_data_std = compute_bands(basis_std, kcoords; n_bands=n_bands, ρ=ρ_ref, tol)
    band_data_mod = compute_bands(basis, kcoords; n_bands=n_bands, ρ=ρ_ref, tol)
    εn_std = n->[εnk[n] for εnk in band_data_std.λ]
    εn_mod = n->[εnk[n] for εnk in band_data_mod.λ]

    # Compute Fermi-Levels
    εF_std = DFTK.compute_occupation(band_data_std.basis, band_data_std.λ).εF
    εF_mod = DFTK.compute_occupation(band_data_mod.basis, band_data_mod.λ).εF
    
    # Return ref_data and mod_data
    (;std_data=(band_data_std, εn_std, εF_std), mod_data=(band_data_mod, εn_mod, εF_mod))
end

function extract_blowup_rate(model)
    (model.model_name=="ModifiedKinetic") &&
        (return model.term_types[1].blowup.p)
    NaN
end
extract_blowup_rate(basis::PlaneWaveBasis) = extract_blowup_rate(basis.model)

function focus_on_band(n, basis_in; ref_data,
                       num_k=100, path_section=ref_data.kpath.paths[1][1:2],
                       tol=1e-4,
                       maxiter=200)
    # Compute zone to zoom on
    kpath = ref_data.kpath
    k_start_label, k_end_label = path_section
    blowup_rate = extract_blowup_rate(basis_in)
    @info "Focusing on band $n between $(k_start_label) and $(k_end_label) with "*
          "$(num_k) points.\n"*"Blow-up rate: $(blowup_rate)"

    # Pre-computations
    kcoords = generate_kpath(kpath.points[k_start_label], kpath.points[k_end_label], num_k)
    ρ_ref = ref_data.scfres.ρ

    # Compute band with higher accuracy between two band diagram points
    # modified kinetic term
    band_data = compute_bands(basis_in, kcoords, n_bands=n, ρ=ρ_ref, tol=tol, maxiter=maxiter)
    εn = [εnk[n] for εnk in band_data.λ]

    # DEBUG
    # subset(tab, n) = [x for (i,x) in enumerate(tab) if rem(i,n)==0]
    # tmp_εn, tmp_kcoords = subset.((εn, kcoords), Ref(debug))
    ∂εn = band_derivative(εn, kcoords)
    ∂2εn = band_derivative(∂εn, kcoords)

    (;path_section, data=(band_data, εn, ∂εn, ∂2εn))
end

# Compute Fermi levels for band diagrams corresponding to the given
# blowup rates.
function compute_fermi_levels(blowup_rates, Ecut::T, n_bands::Int64;
                       ref_data, interval=DefaultInterval,
                       tol=1e-10) where {T<:Real}

    @info "Compute modified kinetic terms bands for low Ecut = $Ecut"
    ref_fft_size = ref_data.scfres.basis.fft_size
    kcoords = ref_data.kcoords
    ρ_ref = ref_data.scfres.ρ

    # Same for each blowup
    εF_list = []
    for blowup_rate in blowup_rates
        blowup = VariableBlowupCHV(blowup_rate; interval)
        basis = ref_data.system.basis(Kinetic(;blowup), Ecut=Ecut, fft_size=ref_fft_size)        
        band_data_mod = compute_bands(basis, kcoords; n_bands=n_bands, ρ=ρ_ref, tol)
        εF = DFTK.compute_occupation(band_data_mod.basis, band_data_mod.λ).εF
        push!(εF_list, (εF, blowup_rate))
    end
    εF_list
end

function test_eigensolver(Ecuts::Vector{T}, blowuprate::T, n_bands::TI;
                          ref_data, interval=DefaultInterval, tol=1e-10,
                          ) where {T<:Real, TI<:Int}
    # Extract ref data
    basis_ref = ref_data.scfres.basis
    fft_size = basis_ref.fft_size
    @assert findmax(Ecuts)[1] ≤ basis_ref.Ecut "Given Ecuts have to be lower that the "*
        "reference Ecut."
    
    # Choose between standard and modified kinetic term.
    blowup = DFTK.BlowupIdentity()
    !isnan(blowuprate) && (blowup = VariableBlowupCHV(blowuprate; interval))

    # Test eigensolver for all given Ecuts    
    res = Dict{T, Any}()
    for Ecut in Ecuts
        @info "Test eigensolver for Ecut=$(Ecut), tol=$(tol) and blow-up rate $(blowuprate)"
        basis = ref_data.system.basis(DFTK.Kinetic(;blowup); Ecut, fft_size)
        H = DFTK.Hamiltonian(basis; ρ = ref_data.scfres.ρ)
        foo = DFTK.diagonalize_all_kblocks(DFTK.lobpcg_hyper, H, n_bands; tol, maxiter=500)
        Bk_size = map(X->size(X,1), foo.X) # Store the number of raws of each Hamiltonian Hk
        res[Ecut] = (; foo.converged, foo.iterations, foo.n_matvec, foo.residual_norms, Bk_size)
    end
    res
end
function test_eigensolver(Ecuts::Vector{T}, blowuprate::T, tols::Vector{T},
                          n_bands::TI;
                          ref_data, file="test_eigensolver",
                          ) where {TI<:Int, T<:Real}
    data = Dict{T, Any}()
    for tol in tols
        data[tol] = test_eigensolver(Ecuts, blowuprate, n_bands; ref_data, tol)
        # Save data in JSON file
        open(io->JSON3.write(io, data, allow_inf=true), file*".json", "w")
    end
    nothing    
end
