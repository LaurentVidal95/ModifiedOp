code_dir = "/home/lvidal/Documents/CERMICS/these/Modified_kinetic_term/code"

const reffile_si = joinpath(code_dir, "../after_modop/bulk_modulus/silicon_PBE/EvsV_SiPBE.json")
const reffile_gr = joinpath(code_dir, "../after_modop/bulk_modulus/graphene_PBE/EvsV_GrPBE.json")

default_curve_params = (; per_of_variation = 5, N_SCF = 51)

# Define usual blowup
interval = [0.85, 0.90] # interpolation interval
blowup_rate = 5/4      # blow-rate of the blow-up function at 1⁻
blowup = VariableBlowupCHV(blowup_rate; interval);

function compute_E₀_vs_V(test_case, Ecut, kgrid, KineticTerm;
                         curve_params=default_curve_params,
                         reference_file)

    @info "Storing data in file $(reference_file)"
    
    # Define a range of lattice constants around equilibrium constant a₀
    system = test_case(; kgrid)
    a₀ = system.a₀
    N_SCF = curve_params.N_SCF
    a_min = (1-curve_params.per_of_variation/100)*a₀
    a_max = (1+curve_params.per_of_variation/100)*a₀
    a_list = LinRange(a_min, a_max, N_SCF)

    # Set SCF params for given test case
    output = Dict{String, Any}()
    blowup_rate = zero(a₀)
    if typeof(KineticTerm.blowup) ≠ BlowupIdentity
        blowup_rate = KineticTerm.blowup.p
    end
    kinlabel = iszero(blowup_rate) ? :standard : :modified
    output["Parameters"] = (; Ecut, kgrid, system.name,
                            KineticTerm=kinlabel,
                            blowup_rate,
                            curve_params...)
    T = eltype(Ecut)
    output["LatticeConstants"] = T[]
    output["Volumes"] = T[]
    output["Energies"] = T[]
    
    # TIME CONSUMING PART
    # Launch SCFs for all constants for given parameters
    progress = Progress(N_SCF, desc="Computing SCFs for $(system.name)...")
    for a in a_list
        scfres = system.scf(; KineticTerm, Ecut, a, n_bands=8,
                            callback=info->nothing)
        V = scfres.basis.model.unit_cell_volume
        E₀ = scfres.energies.total
        push!(output["LatticeConstants"], a)
        push!(output["Volumes"], V)
        push!(output["Energies"], E₀)
        
        # Save result for backup
        open(io->JSON3.write(io, output, allow_inf=true),
             "/home/../tmp/backup.json", "w")
        next!(progress)
    end
    
    # Enter first data in the reference file
    if !isfile(reference_file)
        refdata = Dict{Symbol, Any}()
        refdata[:standard] = Dict{T, Any}()
        refdata[:modified] = Dict{T, Any}()
        refdata[kinlabel][Ecut] = output
        open(io->JSON3.write(io, refdata, allow_inf=true), reference_file, "w")
    else # or add the data to the reference file
        refdata = open(JSON3.read, reference_file)
        newrefdata = copy_nested_dict(refdata)
        newrefdata[kinlabel][Ecut] = output
        open(io->JSON3.write(io, newrefdata, allow_inf=true), reference_file, "w")
    end

    rm("/home/../tmp/backup.json")
    
    output
end

function ∂(tab::Vector{T}, Δ::T) where {T<:Real}
    [(tab[i+1] - tab[i]) / Δ for i in 1:length(tab)-1]
end

# Needed since its hard to add an entry to a JSON dict.
function copy_nested_dict(dict)
    T = eltype(keys(dict))
    new_dict = Dict()
    for key in keys(dict)
        if typeof(dict[key]) <: AbstractDict
            new_dict[key] = copy_nested_dict(dict[key])
        else
            new_dict[key] = dict[key]
        end
    end
    new_dict
end

# # Elastic constant

# """
# Compute the stresses (= 1/Vol dE(a)/da around the a₀ of an obtained SCF solution.
# """
# function compute_stresses_wr_lattice_cste(scfres, a₀::T) where {T<:Real}
#     function HF_energy(a::T) where {T}
#         basis = scfres.basis
#         new_lattice = (a/a₀) .* scfres.basis.model.lattice
#         new_model = Model(basis.model; lattice=new_lattice)
#         new_basis = PlaneWaveBasis(new_model,
#                                    basis.Ecut, basis.fft_size, basis.variational,
#                                    basis.kcoords_global, basis.kweights_global,
#                                    basis.kgrid, basis.kshift, basis.symmetries_respect_rgrid,
#                                    basis.comm_kpts, basis.architecture)
#         ρ = compute_density(new_basis, scfres.ψ, scfres.occupation)
#         energies = energy_hamiltonian(new_basis, scfres.ψ, scfres.occupation;
#                                       ρ, scfres.eigenvalues, scfres.εF).energies
#         energies.total
#     end
#     L = scfres.basis.model.lattice
#     Ω = scfres.basis.model.unit_cell_volume
#     stresses = ForwardDiff.derivative(a -> HF_energy(a), a₀) / Ω
#     # DFTK.symmetrize_stresses(scfres.basis, stresses)
# end

# function compute_stresses_wr_lattice_cste(model::Model, a₀::T; basis_kwargs...) where {T<:Real}
#     basis = PlaneWaveBasis(model; basis_kwargs...)
#     scfres = self_consistent_field(basis; callback=info->nothing)
#     compute_stresses_wr_lattice_cste(scfres, a₀)
# end

# function compute_elastic_constant_wr_lattice_cste(model::Model, a₀::T; basis_kwargs...) where {T<:Real}
#     ForwardDiff.derivative(a -> compute_stresses_wr_lattice_cste(model, a; basis_kwargs...), a₀)
# end
