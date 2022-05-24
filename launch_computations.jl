"""
Main file launching all computations to generate plots.
Estimated time: 5-10 min depending on the systems and parameters.
"""

include("include_utils.jl")

# Choose test case and set parameters
n_bands = 8
blow_up_rate = -1

# system = silicon_PBE(Ecut_ref=Ecut_ref, n_bands=n_bands)
system = graphene_PBE(Ecut_ref=20, kshift=zeros(3), n_bands=n_bands)

# Compute reference band with large Ecut
ref_data = reference_data(system)

# Compute bands for low Ecut with standard and modified kinetic term
# May be a time consuming step
Ecut = 5
plot_data = bandstructure_data(Ecut, n_bands, blow_up_rate; ref_data)

# generate and save plots
!isdir("output") && (mkdir("output"))
savedir="output/$(system.name)"
!isdir(savedir) && (mkdir(savedir))

plot_gm(ha(0.4, blow_up_rate); savedir)      # blow up function
plot_M_EC(ref_data; savedir)   # M_EC standard along kpath
plot_bandstructures(plot_data; ref_data, savedir) # reference, standard and modified terms bandplots

# Focus on a given band
n = 1 # 1 for silicon, 4 for graphene.
band_n_data = focus_on_band(n, blow_up_rate; ref_data, num_k=1000, path_section=["X","U"] )
plot_band(band_n_data; ref_data, savedir)
plot_derivatives_wr_blow_up_rate(n, [-1,-2]; ref_data, num_k=1000, savedir,
                                 path_section=["X","U"]);
