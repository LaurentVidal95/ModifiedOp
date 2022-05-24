include("include_utils.jl")

# Choose test case
Ecut_ref = 10
n_bands = 8
system = silicon_PBE(Ecut_ref=Ecut_ref, n_bands=n_bands)
# graphene_PBE(Ecut_ref=15, kshift=zeros(3))

# Choose blowup function
a = 0.4 # cste to match x^2 in [0,3/4)
ε = -2    # blow_up rate is |⋅|^ε
blow_up_function = y->gm(y, ha(a, ε))

# Compute reference band with large Ecut
ref_data = reference_data(system)

# Compute bands for low Ecut with standard and modified kinetic term
# May be a time consuming step
Ecut = 5
plot_data = compute_band_and_derivatives(Ecut, n_bands, blow_up_function; ref_data)

# generate and save plots
!isdir("output") && (mkdir("output"))
savedir="output/$(system.name)"
!isdir(savedir) && (mkdir(savedir))

plot_gm(ha(a,ε); savedir)      # blow up function
plot_M_EC(ref_data; savedir)   # M_EC standard along kpath
plot_bands(plot_data; ref_data, savedir) # reference, standard and modified terms bandplots

# Focus on a given band
n = 1 # 4 for graphene
plot_band_irregularity(n, plot_data; ref_data, num_k=1000, savedir)
plot_band_derivatives(n, plot_data; ref_data, savedir);
