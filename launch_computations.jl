include("include_utils.jl")

# # Choose test case and set parameters
# Ecut_ref = 10
# n_bands = 8
# blow_up_rate = -2

# system = silicon_PBE(Ecut_ref=Ecut_ref, n_bands=n_bands)
# # graphene_PBE(Ecut_ref=15, kshift=zeros(3))

# # Compute reference band with large Ecut
# ref_data = reference_data(system)

# Compute bands for low Ecut with standard and modified kinetic term
# May be a time consuming step
# Ecut = 5
# plot_data = bandstructure_data(Ecut, n_bands, blow_up_rate; ref_data)

# # generate and save plots
# !isdir("output") && (mkdir("output"))
# savedir="output/$(system.name)"
# !isdir(savedir) && (mkdir(savedir))

# plot_gm(ha(a,ε); savedir)      # blow up function
# plot_M_EC(ref_data; savedir)   # M_EC standard along kpath
# plot_bandstructures(plot_data; ref_data, savedir) # reference, standard and modified terms bandplots

# # Focus on a given band
# n = 1 # 4 for graphene
# plot_band_irregularity(n, plot_data; ref_data, num_k=1000, savedir)
# plot_band_derivatives(n, plot_data; ref_data, savedir);
