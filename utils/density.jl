# ρ = r_to_G(basis, scfres.ρ)
# scatter(vcat(G_energies...), norm.(vcat(test...)), yscale=:log10, ylims=(1e-12,1), xlabel="Energy", ylabel="|ρ|^2")

Hs_norm(basis:: PlaneWaveBasis, ψ_fourier, s) = sqrt( sum((1+norm(G)^2)^s * abs2(ψ_fourier[iG])
                                              for (iG,G) in enumerate(G_vectors_cart(basis))) )

function regularized_density(basis::PlaneWaveBasis, ρ; E_reg=20)
    G_energies = [norm(G)^2/2 for G in G_vectors_cart(basis)][:]
    ρ_reg = zero(ρ)
    ρ_reg[G_energies .≤ E_reg] .= ρ[G_energies .≤ E_reg]
    # Check that the error of approximation is not to high
    println("Error of H1 approximation: $(Hs_norm(basis, ρ_reg .- ρ, 1)*100)%")
    ρ_reg
end
