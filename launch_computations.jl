"""
Main file launching all computations to generate plots.
Estimated time: 5-10 min depending on the systems and parameters.
"""

include("include_utils.jl")

# Choose test case and set parameters
n_bands = 8
Ecut = 5
blow_up_rate = -1
interp_interval = [0.5, 0.75]
system = silicon_PBE(Ecut_ref=15, n_bands=n_bands)

# system = graphene_PBE(Ecut_ref=20, kshift=zeros(3), n_bands=n_bands)

# Compute reference band with large Ecut
# ref_data = reference_data(system)

# Compute bands for low Ecut with standard and modified kinetic term
# May be a time consuming step
bandplot_data = bandstructure_data(Ecut, n_bands, blow_up_rate;
                                   ref_data, interp_interval)

# generate and save plots
# !isdir("../output/other_blow_up_functions") && (mkdir("../output/other_blow_up_functions"))
# savedir="../output/other_blow_up_functions/$(system.name)"
# !isdir(savedir) && (mkdir(savedir))

# plot_gm(blow_up_rate; savedir, interp_interval)      # blow up function
# plot_M_EC(ref_data; savedir)                         # M_EC standard along kpath

# # reference, standard and modified terms bandplots
# plot_bandstructures(bandplot_data; ref_data, savedir) 

# ## Focus on a given band
# n = 1                               # Band to focus on
# path_section=["X", "U"]             # Section of the k_path to focus on
# num_k=1000                          # Resolution of the band
# tol=1e-5                            # Tolerance for the eigensolver for the Hk
# blow_up_rates = [-1, -2]      # All blow-up rates to test

# # Reference fft grid so that every one has the same reference ground state density.
# ref_fft_size = ref_data.scfres.basis.fft_size

# basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
# basis_mod(blow_up_rate) = ref_data.system.basis(
#     ModifiedKinetic(blow_up=y->gm(y, ha(0.4, blow_up_rate))),
#     Ecut=Ecut, fft_size=ref_fft_size)

# # Computation of the bands
# band_ref_data = focus_on_band(n, ref_data.scfres.basis; ref_data,
#                                 path_section, num_k, tol)
# band_std_data = focus_on_band(n, basis_std; ref_data,
#                               path_section, num_k, tol)
# bands_mod_data = focus_on_band.(Ref(n), basis_mod.(blow_up_rates);
#                                   ref_data, path_section, num_k, tol)

# # Plot
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=0)
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=1)
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=2)
