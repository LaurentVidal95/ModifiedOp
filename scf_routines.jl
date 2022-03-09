# TODO
# Plot derivative with finite diff

using Pkg; Pkg.activate("/home/lvidal/programs/DFTK_perso.jl"); using DFTK;
using Plots
using Unitful
using StaticArrays

include("g_function.jl")
include("k_path_perso.jl")
include("energy_bands.jl")

### Terms
reference_PBE_terms = [Kinetic(), AtomicLocal(), AtomicNonlocal(),
                      Ewald(), PspCorrection(), Hartree(), Xc([:gga_x_pbe, :gga_c_pbe])]
modified_PBE_terms(;g=y->y'y) = [RegularizedKinetic(;g=g),
                           AtomicLocal(), AtomicNonlocal(), Ewald(),
                           PspCorrection(), Hartree(), Xc([:gga_x_pbe, :gga_c_pbe])]

"""
    Graphene
"""
function model_PBE_graphene(terms)
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
    Model(lattice; atoms=atoms, terms=terms)
end

function scf_graphene(;n_bands=13, terms=reference_PBE_terms,
                      Ecut=15, kgrid=[8,8,1], kshift=zeros(3))
    @show n_bands, E_cut, kgrid, kshift
    model = model_PBE_graphene(terms)
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid, kshift=kshift)
    scfres = self_consistent_field(basis, n_bands=n_bands);
end
