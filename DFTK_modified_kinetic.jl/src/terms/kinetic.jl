"""
Kinetic energy: 1/2 sum_n f_n ∫ |∇ψn|^2.
"""
struct Kinetic
    scaling_factor::Real
end
Kinetic(; scaling_factor=1) = Kinetic(scaling_factor)
(kin::Kinetic)(basis) = TermKinetic(basis, kin.scaling_factor)
function Base.show(io::IO, kin::Kinetic)
    fac = isone(kin.scaling_factor) ? "" : ", scaling_factor=$scaling_factor"
    print(io, "Kinetic($fac)")
end

struct TermKinetic <: Term
    scaling_factor::Real  # scaling factor, absorbed into kinetic_energies
    kinetic_energies::Vector{<:AbstractVector}  # kinetic energy 1/2|G+k|^2 for every kpoint
end
function TermKinetic(basis::PlaneWaveBasis{T}, scaling_factor) where {T}
    kinetic_energies = [[T(scaling_factor) * sum(abs2, Gk) / 2
                         for Gk in Gplusk_vectors_cart(basis, kpt)]
                        for kpt in basis.kpoints]
    TermKinetic(T(scaling_factor), kinetic_energies)
end

@timing "ene_ops: kinetic" function ene_ops(term::TermKinetic, basis::PlaneWaveBasis{T},
                                            ψ, occ; kwargs...) where {T}
    ops = [FourierMultiplication(basis, kpoint, term.kinetic_energies[ik])
           for (ik, kpoint) in enumerate(basis.kpoints)]
    isnothing(ψ) && return (E=T(Inf), ops=ops)

    E = zero(T)
    for (ik, k) in enumerate(basis.kpoints)
        for iband = 1:size(ψ[1], 2)
            ψnk = @views ψ[ik][:, iband]
            E += (basis.kweights[ik] * occ[ik][iband]
                  * real(dot(ψnk, Diagonal(term.kinetic_energies[ik]), ψnk)))
        end
    end
    E = mpi_sum(E, basis.comm_kpts)

    (E=E, ops=ops)
end

struct ModifiedKinetic
    scaling_factor::Real
    blow_up_function # gm function
end
# Defaut is standard Kinetic term but with strict inequality in the choice of G_vectors
ModifiedKinetic(; scaling_factor=1, blow_up=x->x'x) = ModifiedKinetic(scaling_factor, blow_up)
(kin::ModifiedKinetic)(basis) = TermKinetic(basis, kin.scaling_factor, kin.blow_up_function)

function TermKinetic(basis::PlaneWaveBasis{T}, scaling_factor, gm) where {T}
    Ecut = basis.Ecut
    Gplusk_strict_cart(basis, kpt) = [Gpk for Gpk in Gplusk_vectors_cart(basis, kpt)
                                      if (norm(Gpk)^2/2 < Ecut)]
    kinetic_energies = [[T(scaling_factor) * Ecut*gm(norm(Gk)/√(2*Ecut))
                         for Gk in Gplusk_strict_cart(basis, kpt)]
                        for kpt in basis.kpoints]
    TermKinetic(T(scaling_factor), kinetic_energies)
end