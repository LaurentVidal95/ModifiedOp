# Elastic constant

"""
Compute the stresses (= 1/Vol dE(a)/da around the a₀ of an obtained SCF solution.
"""
function compute_stresses_wr_lattice_cste(scfres, a₀::T) where {T<:Real}
    function HF_energy(a::T) where {T}
        basis = scfres.basis
        new_lattice = (a/a₀) .* scfres.basis.model.lattice
        new_model = Model(basis.model; lattice=new_lattice)
        new_basis = PlaneWaveBasis(new_model,
                                   basis.Ecut, basis.fft_size, basis.variational,
                                   basis.kcoords_global, basis.kweights_global,
                                   basis.kgrid, basis.kshift, basis.symmetries_respect_rgrid,
                                   basis.comm_kpts, basis.architecture)
        ρ = compute_density(new_basis, scfres.ψ, scfres.occupation)
        energies = energy_hamiltonian(new_basis, scfres.ψ, scfres.occupation;
                                      ρ, scfres.eigenvalues, scfres.εF).energies
        energies.total
    end
    L = scfres.basis.model.lattice
    Ω = scfres.basis.model.unit_cell_volume
    stresses = ForwardDiff.derivative(a -> HF_energy(a), a₀) / Ω
    # DFTK.symmetrize_stresses(scfres.basis, stresses)
end

function compute_stresses_wr_lattice_cste(model::Model, a₀::T; basis_kwargs...) where {T<:Real}
    basis = PlaneWaveBasis(model; basis_kwargs...)
    scfres = self_consistent_field(basis; callback=info->nothing)
    compute_stresses_wr_lattice_cste(scfres, a₀)
end

function compute_elastic_constant_wr_lattice_cste(model::Model, a₀::T; basis_kwargs...) where {T<:Real}
    ForwardDiff.derivative(a -> compute_stresses_wr_lattice_cste(model, a; basis_kwargs...), a₀)
end
