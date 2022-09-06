"""
Main file launching all computations to generate plots.
Estimated time: 5-20 min depending on the systems and parameters.
The parameters influencing time are:
    - `bandplot_res` which defines the resolution of kpoint between two
      symmetry points of the band diagram (Beware: bandplot_res is not the total number
      of k-points).
    - `single_band_res` which is the number of k-points used to compute a single
      band and its finite difference derivatives. It should be high enough to
      avoid oscillations due to the finite difference approximation, especially
      at order 2.

To launch computation simply call launch_computation, e.g. for silicon with extra
rough parameters (computation takes 1min max) define an output directory "../silicon" and
call in terminal:

``
launch_computations(silicon, blowup; bandplot_res=50, single_band_res=100, output_dir="../silicon",
Ecut, n_bands)
``

All functions are defined in the `utils` directory.
"""

include("include_utils.jl")

# Set general parameters, OK for both systems
n_bands = 8
Ecut = 5

# Define blow-up
interval = [0.85, 0.90] # interpolation interval
blowup_rate = 3//2
blowup = VariableBlowupCHV(blowup_rate; interval)

# Define system
silicon = silicon_PBE(; Ecut_ref=20, n_bands, kgrid=[10,10,10])
graphene = graphene_PBE(; Ecut_ref=20, kshift=zeros(3), n_bands, kgrid=[12,12,1])

# Now only launch in terminal the following function
function launch_computations(system, blowup; bandplot_res=200, single_band_res=2000,
                             output_dir="", Ecut, n_bands)
    # Compute reference band with large Ecut
    @info "Computing reference ground state density"
    ref_data = reference_data(system; k_path_res=bandplot_res)

    # Compute bands for low Ecut with standard and modified kinetic term
    # May be a time consuming step
    bandplot_data = bandstructure_data(Ecut, n_bands, blowup; ref_data, interval)

    # generate output directories
    data_dir = isempty(output_dir) ? "" : joinpath(output_dir, "data")
    (!isempty(data_dir) && !isdir(data_dir)) && (mkdir(data_dir))

    plot_dir = isempty(output_dir) ? "" : joinpath(output_dir, "plots")
    (!isempty(plot_dir) && !isdir(plot_dir)) && (mkdir(plot_dir))

    # Plots on whole band diagram
    plot_M_EC(ref_data; plot_dir)
    p = plot_blowup_function(blowup; plot_dir)
    plot_bandstructures(bandplot_data; ref_data, plot_dir)

    # Focus on a given band
    n = system.n                          # Band to focus on
    path_section = system.path_section    # Section of the k_path to focus on
    num_k = single_band_res               # Resolution of the band
    tol = 1e-12                           # Tolerance for the eigensolver for the Hk
    maxiter = 200                         # max_iterations for lobpcg solver. default=100.
    blow_up_rates = [1/2, 3/2, 5/2]       # All blow-up rates to test

    # Reference fft grid so that every one has the same reference ground state density.
    ref_fft_size = ref_data.scfres.basis.fft_size

    basis_std = ref_data.system.basis(Kinetic(), Ecut=Ecut, fft_size=ref_fft_size)
    basis_mod(blowup_rate) = ref_data.system.basis(
        Kinetic(; blowup=VariableBlowupCHV(blowup_rate; interval)),
        Ecut=Ecut, fft_size=ref_fft_size)

    # Computation of the bands
    band_ref_data = focus_on_band(n, ref_data.scfres.basis; ref_data,
                                  path_section, num_k, tol)
    band_std_data = focus_on_band(n, basis_std; ref_data,
                                  path_section, num_k, tol)
    bands_mod_data = focus_on_band.(Ref(n), basis_mod.(blow_up_rates);
                                    ref_data, path_section, num_k, tol, maxiter)

    # save raw data
    if !(isempty(data_dir))
        save_band_data(band_ref_data, joinpath(data_dir, "ref_data"))
        save_band_data(band_std_data, joinpath(data_dir, "std_data"))
        for i in 1:length(blow_up_rates)
            save_band_data(bands_mod_data[i], joinpath(data_dir,
                                                       "mod_data_$(blow_up_rates[i])"))
        end
    end

    # Plot
    plot_band(band_ref_data, band_std_data, bands_mod_data;
              ref_data, plot_dir, i_derivative=0)
    plot_band(band_ref_data, band_std_data, bands_mod_data;
              ref_data, plot_dir, i_derivative=1)
    plot_band(band_ref_data, band_std_data, bands_mod_data;
              ref_data, plot_dir, i_derivative=2)
end

## EXAMPLES: launch computation for silicon
## bandplot_res and single_band_res have been set to avoid long computation time.
## The results of the paper are displayed for bandplot_res=300 and single_band_res=4000.

# mkdir("../silicon_PBE")
# launch_computations(silicon, blowup; bandplot_res=100, single_band_res=100, output_dir="", Ecut, n_bands)
