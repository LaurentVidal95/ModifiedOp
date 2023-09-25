function FD_derivative(tab, Δ, deg=1)
    @assert (deg ≥ 0)
    (deg==0) && return tab
    # centered FD for O(h^2) error
    ∂tab = [(tab[i+1] - tab[i-1]) / (2*Δ) for i in 2:length(tab)-1]
    FD_derivative(∂tab, Δ, deg-1)
end

function HF_energy_derivative(a, system, Ecut, KineticTerm;
                              hacklist=[])
    scfres = system.scf(; Ecut, a, tol=1e-10, n_bands=8, callback=identity)

    # hack to extract all information while only calling second derivative
    # routine. Allow to compute E, dE and d^2E in a single call.
    push!(hacklist, filter_dual(a))
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
        # Extract total energy value from loop
        E = energies.total
        push!(hacklist, filter_dual(E))
        E
    end
    dE = ForwardDiff.derivative(HF_energy, a)
    push!(hacklist, filter_dual(dE))

    dE
end
function HF_energy_second_derivative(a₀, system, Ecut, KineticTerm;
                                     return_energy_and_derivative=false)
    hacklist = []
    d2E = ForwardDiff.derivative(
        a->HF_energy_derivative(a, system, Ecut, KineticTerm; hacklist), a₀
    )
    push!(hacklist, d2E)
    (return_energy_and_derivative) && (return hacklist)
    d2E
end

# TMP instruction. See TODO 3)
# progress = Progress(length(a_list), desc="Computing HF dev...")
# output = Dict{String, Any}()
# output["LatticeConstants"] = a_list
# output["Energies"] = []
# output["FirstDerivative"] = []
# output["SecondDerivative"] = []
# output["Parameters"] = (;Ecut=30, blowuprate=0.0, system="Silicon_PBE", kgrid=[5,5,5])
# for a in a_list
#     _, E, dE, d2E = HF_energy_second_derivative(a, Si, 30, Kinetic();
#                                                 return_energy_and_derivative=true)
#     push!(output["Energies"], E)
#     push!(output["FirstDerivative"], dE)
#     push!(output["SecondDerivative"], d2E)
#     open(io->JSON3.write(io, output, allow_inf=true), "HF_dev_backup.json", "w")
#     next!(progress)
# end
