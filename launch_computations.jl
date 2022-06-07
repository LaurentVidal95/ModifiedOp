"""
Main file launching all computations to generate plots.
Estimated time: 5-10 min depending on the systems and parameters.
"""

include("include_utils.jl")

# Choose test case and set parameters
n_bands = 8

# system = silicon_PBE(Ecut_ref=15, n_bands=n_bands)
# # system = graphene_PBE(Ecut_ref=20, kshift=zeros(3), n_bands=n_bands)

# # Compute reference band with large Ecut
# ref_data = reference_data(system)

# # Compute bands for low Ecut with standard and modified kinetic term
# # May be a time consuming step
Ecut = 5
# bandplot_data = bandstructure_data(Ecut, n_bands, blow_up_rate; ref_data)

# generate and save plots
!isdir("../output/confirm_proof") && (mkdir("../output/confirm_proof"))
savedir="../output/confirm_proof/$(system.name)"
!isdir(savedir) && (mkdir(savedir))

# plot_gm(ha(0.4, blow_up_rate); savedir)      # blow up function
# plot_M_EC(ref_data; savedir)   # M_EC standard along kpath
# # reference, standard and modified terms bandplots
# plot_bandstructures(bandplot_data; ref_data, savedir) 

# Focus on a given band
n = 1
num_k=100
blow_up_rates = [-1//2, -5//2]
path_section=["X", "U"]
ref_fft_size = ref_data.scfres.basis.fft_size

# Compute std and mod basis on same fft_grid as reference to be able to take
# same GS density ρ_ref as reference.
basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
basis_mod(blow_up_rate) = ref_data.system.basis(
    ModifiedKinetic(blow_up=y->gm(y, ha(0.4, blow_up_rate))),
    Ecut=Ecut, fft_size=ref_fft_size)

# Reference band
# band_ref_data = focus_on_band(n, ref_data.scfres.basis; ref_data,
#                                 path_section, num_k)
# band_std_data = focus_on_band(n, basis_std; ref_data, path_section)
# bands_mod_data = focus_on_band.(Ref(n), basis_mod.(blow_up_rates);
#                                   ref_data, path_section)

plot_band(band_ref_data, band_std_data, bands_mod_data;
          ref_data, savedir, i_derivative=2)
# band_n_data = focus_on_band(n, blow_up_rate; ref_data, num_k=1000, path_section=["X","U"] )
# plot_band(band_n_data; ref_data, savedir)
# plot_derivatives_wr_blow_up_rate(n, [-0.5, -1.5,-2.5]; ref_data, num_k=1000, savedir,
#                                  path_section=["X","U"]);
