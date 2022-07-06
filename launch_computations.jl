"""
Main file launching all computations to generate plots.
Estimated time: 5-20 min depending on the systems and parameters.
"""

include("include_utils.jl")

# Choose test case and set parameters
n_bands = 8
Ecut = 5
# Blow_up_function
interval = [0.97, 0.98]
blow_up_rate = -3//2
blow_up_function = gm(blow_up_rate; interval)

# System
system = silicon_PBE(Ecut_ref=15, n_bands=n_bands)
# system = graphene_PBE(Ecut_ref=20, kshift=zeros(3), n_bands=n_bands)

# Compute reference band with large Ecut
ref_data = reference_data(system)

# Compute bands for low Ecut with standard and modified kinetic term
# May be a time consuming step
# bandplot_data = bandstructure_data(Ecut, n_bands, blow_up_function;
#                                    ref_data, interval)

# # generate and save plots
!isdir("../output/confirm_proof") && (mkdir("../output/confirm_proof"))
savedir="../output/confirm_proof/$(system.name)"
!isdir(savedir) && (mkdir(savedir))

# Plots on whole band diagram
plot_M_EC(ref_data; savedir)
plot_blow_up_function(blow_up_function; savedir)
# plot_bandstructures(bandplot_data; ref_data, savedir) 


# # Focus on a given band
n=1  #n = 2                         # Band to focus on
path_section=["X", "U"]             # Section of the k_path to focus on
num_k=2000                           # Resolution of the band
tol=1e-12                           # Tolerance for the eigensolver for the Hk
maxiter=200                         # max_iterations for lobpcg solver. default=100.
blow_up_rates = [-1/2, -3/2, -5/2]  # All blow-up rates to test

# Reference fft grid so that every one has the same reference ground state density.
ref_fft_size = ref_data.scfres.basis.fft_size

basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
basis_mod(blow_up_rate) = ref_data.system.basis(
    ModifiedKinetic(blow_up=gm(blow_up_rate; interval)),
    Ecut=Ecut, fft_size=ref_fft_size)

# Computation of the bands
band_ref_data = focus_on_band(n, ref_data.scfres.basis; ref_data,
                                path_section, num_k, tol)
band_std_data = focus_on_band(n, basis_std; ref_data,
                              path_section, num_k, tol)
bands_mod_data = focus_on_band.(Ref(n), basis_mod.(blow_up_rates);
                                  ref_data, path_section, num_k, tol, maxiter)

# Plot
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=0)
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=1)
# plot_band(band_ref_data, band_std_data, bands_mod_data;
#           ref_data, savedir, i_derivative=2)

"""
Silicon
dev 1:
lens!([x_axis[437],x_axis[442]], [-0.0002, 0.00],
                inset = (1, bbox(0.5, 0.5, 0.4, 0.4)), subplot=2)
xticks!(p.subplots[2], :none)
yticks!(p.subplots[2], :none)

dev 2:
 lens!([x_axis[435],x_axis[442]], [-0.00015, -0.0001], 
                inset = (1, bbox(0.5, 0.5, 0.4, 0.4)), subplot=2)
xticks!(p.subplots[2], :none)
yticks!(p.subplots[2], :none)
"""
