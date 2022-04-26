using Pkg; Pkg.activate("DFTK_modified_kinetic.jl"); using DFTK;
using Plots, Measures
using Unitful, UnitfulAtomic
using DelimitedFiles
using LaTeXStrings

# Include utils dir
include.(joinpath.(Ref("utils"), readdir("utils/")))

"""
Terms
"""
PBE_terms(KineticTerm) = [KineticTerm, AtomicLocal(), AtomicNonlocal(),
                       Ewald(), PspCorrection(), Hartree(), Xc([:gga_x_pbe, :gga_c_pbe])]
# LDA_terms(KinteticTerm) = [KineticTerm, AtomicLocal(), AtomicNonLocal(),
#                        Ewald(), PspCorrection(), Hartree(), Xc([:lda_xc_teter93])]

"""
Graphene
"""
function graphene_PBE(; Ecut_ref=15, basis_kwargs...)
    # Defines graphene structure and hamiltonian terms
    function model_PBE_graphene(KineticTerm;
                     a=austrip(1u"Å")*2.641) # Lattice constant of graphene in bohr
        # Define lattice
        a_1 = [a; 0; 0];
        rot_120_deg = [-1/2   -√3/2   0;
                       √3/2    -1/2   0;
                       0         0    0]
        a_2 = rot_120_deg*a_1; a_3 = [0; 0; 20]
        lattice = hcat(a_1, a_2, a_3)
        # Define elements
        C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
        atoms = [C => [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]]
        # PBE functional
        model_name="custom"
        (KineticTerm isa RegularizedKinetic) && (model_name="RegularizedKinetic")
        Model(lattice; atoms=atoms, terms=PBE_terms(KineticTerm), model_name=model_name)
    end

    # Construct a plane-wave basis given a kinetic term using model_PBE_graphene
    function basis_PBE_graphene(KineticTerm; Ecut=Ecut_ref, basis_kwargs...)
        model = model_PBE_graphene(KineticTerm)
        PlaneWaveBasis(model; Ecut=Ecut, basis_kwargs...)
    end

    # Scf using the above functions.Sets the defaut Ecut, number of kpoints and bands.
    function scf_graphene(;n_bands=n_bands, KineticTerm=Kinetic(), Ecut=Ecut_ref)
        basis = basis_PBE_graphene(KineticTerm; Ecut=Ecut, basis_kwargs...)
        self_consistent_field(basis, n_bands=n_bands);
    end
    (;scf=scf_graphene, basis=basis_PBE_graphene, model=model_PBE_graphene)
end

"""
Silicon
"""
function silicon_PBE(; Ecut_ref=15, basis_kwargs)
    # Defines graphene structure and hamiltonian terms
    function model_PBE_silicon(KineticTerm;
                               a=10.26) # Lattice constant of silicon in bohr
        lattice = a / 2 * [[0 1 1.];
                           [1 0 1.];
                           [1 1 0.]]
        Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
        atoms = [Si => [ones(3)/8, -ones(3)/8]]
        # PBE functional
        model_name="custom"
        (KineticTerm isa RegularizedKinetic) && (model_name="RegularizedKinetic")
        Model(lattice; atoms=atoms, terms=PBE_terms(KineticTerm), model_name=model_name)
    end

    # Construct a plane-wave basis given a kinetic term using model_PBE_graphene
    function basis_PBE_silicon(KineticTerm; Ecut=Ecut_ref, basis_kwargs...)
        model = model_PBE_silicon(KineticTerm)
        PlaneWaveBasis(model; Ecut=Ecut, basis_kwargs...)
    end

    # Scf using the above functions.Sets the defaut Ecut, number of kpoints and bands.
    function scf_silicon(;n_bands=n_bands, KineticTerm=Kinetic(), Ecut=Ecut_ref)
        basis = basis_PBE_silicon(KineticTerm; Ecut=Ecut, basis_kwargs...)
        self_consistent_field(basis, n_bands=n_bands);
    end
    (;scf=scf_graphene, basis=basis_PBE_silicon, model=model_PBE_silicon)
end

# """
# Silicon cohen-bergstresser
# """
# function model_cohen_bergstresser_silicon(KineticTerm) 
#     Si = ElementCohenBergstresser(:Si)
#     lattice = Si.lattice_constant / 2 .* [[0 1 1.]; [1 0 1.]; [1 1 0.]]
#     atoms = [Si => [ones(3)/8, -ones(3)/8]];
#     # Change model name
#     model_name="custom";
#     (KineticTerm isa RegularizedKinetic) && (model_name="RegularizedKinetic")
#     model = Model(lattice; atoms=atoms, terms=[KineticTerm, AtomicLocal()], model_name=model_name)
# end

# function solve_cohen_bergstresser_silicon(; KineticTerm=Kinetic(),
#                                           n_bands=6,
#                                           Ecut=10.0, kgrid=[2,2,2], kshift=zeros(3))
#     model = model_cohen_bergstresser_silicon(KineticTerm)
#     basis = PlaneWaveBasis(model, Ecut=Ecut, kgrid=kgrid, kshift=kshift)
#     # solve linear problem
#     ham = Hamiltonian(basis)
#     eigres=DFTK.diagonalize_all_kblocks(DFTK.lobpcg_hyper, ham, n_bands)

#     # return as a pseudo_scfres
#     (ψ=eigres.X, eigenvalues=eigres.λ, basis=basis, diagonalization=eigres)
# end
