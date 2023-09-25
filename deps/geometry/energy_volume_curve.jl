code_dir = "/home/lvidal/Documents/CERMICS/these/Modified_kinetic_term/code"

const reffile_si = joinpath(code_dir, "../after_modop/data/EvsV_SiPBE.json")
const reffile_gr = joinpath(code_dir, "../after_modop/data/EvsV_GrPBE.json")

default_curve_params = (; per_of_variation = 5, N_SCF = 51)

# Define usual blowup
interval = [0.85, 0.90] # interpolation interval
blowup_rate = 5/4      # blow-rate of the blow-up function at 1⁻
blowup = VariableBlowupCHV(blowup_rate; interval);

function compute_E₀_vs_lattice_cste(test_case, Ecut, kgrid, KineticTerm;
                                    curve_params=default_curve_params,
                                    backupdir="/tmp",
                                    reference_file)
    system = test_case()
    a₀ = system.a₀
    N_SCF = curve_params.N_SCF
    a_min = (1-curve_params.per_of_variation/100)*a₀
    a_max = (1+curve_params.per_of_variation/100)*a₀
    a_list = LinRange(a_min, a_max, N_SCF)

    compute_E₀_vs_lattice_cste(test_case, Ecut, kgrid, a_list, KineticTerm;
                               reference_file, backupdir)
end

function compute_E₀_vs_lattice_cste(test_case, Ecut, kgrid, a_list, KineticTerm;
                                    reference_file,
                                    backupdir="/tmp")

    # Define a range of lattice constants around equilibrium constant a₀
    system = test_case(; kgrid)
    a₀ = system.a₀
    N_SCF = length(a_list)

    # Set SCF params for given test case
    output = Dict{String, Any}()
    blowup_rate = zero(a₀)
    if typeof(KineticTerm.blowup) ≠ BlowupIdentity
        blowup_rate = KineticTerm.blowup.p
    end
    kinlabel = iszero(blowup_rate) ? :standard : :modified
    output["Parameters"] = (; Ecut, kgrid, system.name,
                            KineticTerm=kinlabel,
                            blowup_rate
                            )
    T = eltype(Ecut)
    output["LatticeConstants"] = T[]
    output["Volumes"] = T[]
    output["Energies"] = T[]

    # TIME CONSUMING PART
    # Launch SCFs for all constants for given parameters
    @info "Compute E=f(a) curve for Ecut=$(Ecut) and blowup rate $(blowup_rate)"
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
             joinpath(backupdir, "backup.json"), "w")
        next!(progress)
    end

    # Store output data and remove backup
    store_energy_curve_data(output, reference_file)
    rm(joinpath(backupdir, "backup.json"))

    output
end

function read_E_vs_lattice_constant_curve(datafile, blowup_rate::T, Ecut::T) where T
    data = open(JSON3.read, datafile)
    a_list = data[blowup_rate][Ecut][:LatticeConstants]
    E₀_list = data[blowup_rate][Ecut][:Energies]
    a_list, E₀_list
end

### Find discontinuity on the curve.
function compute_PW_lengths(system, lattice_constants, Ecut)
    map(lattice_constants) do a
        basis = system.basis(Kinetic(); Ecut, a)
        length(G_vectors(basis))
    end
end

function spot_PW_length_jump(test_case, targeted_precision, Ecut;
                             N_lattice_constants=50,
                             per_of_variation=10)
    # First rough
    system = test_case(;)
    a₀ = system.a₀

    a_min = (1-per_of_variation/100)*a₀
    a_max = (1+per_of_variation/100)*a₀
    a_list = LinRange(a_min, a_max, N_lattice_constants)
    Δa = a_list[2] - a_list[1]

    maxiter = 10
    iter = zero(maxiter)
    i_jump = nothing
    while (Δa > targeted_precision) && (iter ≤ maxiter)
        PW_lengths = compute_PW_lengths(system, a_list, Ecut)
        i_jump = findfirst(i->!iszero(PW_lengths[i+1] - PW_lengths[i]),
                           1:length(PW_lengths)-1)
        isnothing(i_jump) && error("No jump detected")
        a_list = LinRange(a_list[i_jump], a_list[i_jump+1], N_lattice_constants)
        Δa = a_list[2] - a_list[1]
        iter += 1
    end
    (iter > maxiter) && @warn "Maximum iteration reached"
    PW_lengths = compute_PW_lengths(system, a_list, Ecut)
    i_jump = findfirst(i->!iszero(PW_lengths[i+1] - PW_lengths[i]),
                       1:length(PW_lengths)-1)
    a_list[i_jump]
end
