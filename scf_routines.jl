using Pkg; Pkg.activate("/home/lvidal/programs/DFTK_perso.jl"); using DFTK;
using Plots
using Unitful
using StaticArrays

# Include utils dir
include.(joinpath.(Ref("utils"), readdir("utils/")))

"""
Terms
"""
PBE_terms(KineticTerm) = [KineticTerm, AtomicLocal(), AtomicNonlocal(),
                          Ewald(), PspCorrection(), Hartree(), Xc([:gga_x_pbe, :gga_c_pbe])]

"""
Graphene
"""
function model_PBE_graphene(KineticTerm)
    ang_to_bohr = 1.88973; a = ang_to_bohr*2.641; # from wiki
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

function basis_PBE_graphene(KineticTerm; Ecut=15, kgrid=[8,8,1], kshift=zeros(3))
    model = model_PBE_graphene(KineticTerm)
    PlaneWaveBasis(model, Ecut=Ecut, kgrid=kgrid, kshift=kshift)
end

function scf_graphene(;n_bands=13, KineticTerm=Kinetic(),
                      Ecut=15, kgrid=[8,8,1], kshift=zeros(3))
    @show n_bands, Ecut, kgrid, kshift
    basis = basis_PBE_graphene(KineticTerm, Ecut=Ecut, kgrid=kgrid, kshift=kshift)
    self_consistent_field(basis, n_bands=n_bands);
end

"""
Silicon cohen-bergstresser
"""
function model_cohen_bergstresser_silicon(KineticTerm) 
    Si = ElementCohenBergstresser(:Si)
    lattice = Si.lattice_constant / 2 .* [[0 1 1.]; [1 0 1.]; [1 1 0.]]
    atoms = [Si => [ones(3)/8, -ones(3)/8]];
    # Change model name
    model_name="custom";
    (KineticTerm isa RegularizedKinetic) && (model_name="RegularizedKinetic")
    model = Model(lattice; atoms=atoms, terms=[KineticTerm, AtomicLocal()], model_name=model_name)
end

function solve_cohen_bergstresser_silicon(; KineticTerm=Kinetic(),
                                          n_bands=6,
                                          Ecut=10.0, kgrid=[2,2,2], kshift=zeros(3))
    model = model_cohen_bergstresser_silicon(KineticTerm)
    basis = PlaneWaveBasis(model, Ecut=Ecut, kgrid=kgrid, kshift=kshift)
    # solve linear problem
    ham = Hamiltonian(basis)
    eigres=DFTK.diagonalize_all_kblocks(DFTK.lobpcg_hyper, ham, n_bands)

    # return as a pseudo_scfres
    (ψ=eigres.X, eigenvalues=eigres.λ, basis=basis, diagonalization=eigres)
end
