# Default plot parameters
pyplot()
myblue=RGB(95/255,133/255,255/255)
mygreen=RGB(60/255,195/255,26/255)
default(fontfamily="serif",
        linewidth=0.75, framestyle=:box,
        label=nothing, grid=:true,
        linecolor=myblue,
        tickfontsize=12,
        titlefontsize=15,
        legendfontsize=12,
        guidefontsize=15)

"""
Plot gm blow up function on [0,1) given third part ha.
"""
function plot_gm(blow_up_rate; interp_interval=[0.7, 0.75], savedir="")    
    x_axis = LinRange(0, 0.99, 1000)

    # Plot
    p = plot(x_axis, abs.(gm.(x_axis, optimal_ha(blow_up_rate; interp_interval);
                              interp_interval)),
             label=:none, xlims=[0,1])
    
    vline!(interp_interval, label=L"C^2"*" spline interpolation",
           linestyle=:dash, linecolor=mygreen)
    vline!([1], label=:none)
    plot!(x->x^2, label=L"x\mapsto x^2", linestyle=:dash, linecolor=:black)

    # Plot parameters
    plot!(legendfontsize=12, legend=:topleft)
    plot!(size=(750,500))
    xlabel!(L"x")
    ylabel!(L"\mathscr{G}(x)")
    title!("Blow up rate: "*latexstring("|\\cdot|^{$(blow_up_rate)}"))
    ylims!(0,5)

    # Save plot if needed
    !isempty(savedir) && (savefig(p, joinpath(savedir,"blow_up_function.pdf")))

    p
end

"""
Plot M_EC along reference k-path
"""
function plot_M_EC(ref_data; savedir="")
    basis = ref_data.band_data.basis
    kpath = ref_data.kpath
    tmp_plot = DFTK.plot_band_data(ref_data.band_data; ref_data.scfres.εF, kpath.klabels)
    ticks = only(xticks(tmp_plot))

    function compute_M_EC_along_path(basis)
        M_EC = [length(G_vectors(basis, kpt)) for kpt in basis.kpoints]
        mean = round(Int64, sum(M_EC) / length(M_EC))
        num_k = length(M_EC)
        jumps = [abs(M_EC[i+1] - M_EC[i]) for i in 1:num_k-1]
        (M_EC .- mean), mean, jumps, findmax(jumps)
    end
    dist_to_mean, mean, _, _ = compute_M_EC_along_path(basis)

    # plot
    x_axis = LinRange(ticks[1][1], ticks[1][end], length(dist_to_mean))
    p = plot(x_axis, dist_to_mean, linewidth=1.5, label=:none)
             # markersize=6, markershape=:cross,

    title!(split(ref_data.system.name,"_",keepempty=false)[1])
    xlabel!(L"\mathrm{wave\; vector\;}k ")
    ylabel!(L"M_{E_c}(k)-\mathrm{mean}(M_{E_c})")
    plot!(size=(800,500))
    xticks!(ticks...)

    !isempty(savedir) && (savefig(p, joinpath(savedir, "M_EC.pdf")))

    p
end

"""
Plot bandstructures with respective standard and modified kinetic terms
"""
function plot_bandstructures(plot_data; ref_data, savedir="")
    # Extract data
    ρ_ref = ref_data.scfres.ρ
    band_data_ref = ref_data.band_data
    band_data_std = plot_data.std_data[1]
    band_data_mod = plot_data.mod_data[1]

    kpath = ref_data.kpath
    kcoords = kpath.kcoords
    num_k = length(kcoords)

    function define_ylims(data, εF)
        y_min = findmin(data.λ)[1][1]
        y_max = findmax([findmax(λi)[1] for λi in data.λ])[1]
        (y_min, y_max) .- εF
    end

    # Reference
    εF = ref_data.scfres.εF
    p_ref = DFTK.plot_band_data(band_data_ref; εF, kpath.klabels)
    ylims!(p_ref, define_ylims(band_data_ref, εF))
    plot!(p_ref, size=(800,500))
    ylabel!(p_ref, L"\varepsilon_{n,k}-\varepsilon_f\;(\mathrm{hartree})")
    xlabel!(p_ref," ")
    title!(p_ref,"Reference band diagram")

    # Standard
    εF = DFTK.fermi_level(band_data_std.basis, band_data_std.λ)
    p_std = DFTK.plot_band_data(band_data_std; εF=εF, kpath.klabels)
    plot!(p_std, size=(800,500))
    ylims!(p_std, define_ylims(band_data_std, εF))
    xlabel!(p_std, " ")
    ylabel!(p_std, L"\varepsilon_{n,k}^{E_c}-\varepsilon_f\;(\mathrm{hartree})")
    title!(p_std,L"\mathrm{Standard\; kinetic\; term}")

    # Modified
    εF = DFTK.fermi_level(band_data_mod.basis, band_data_mod.λ)
    p_mod = DFTK.plot_band_data(band_data_mod; εF, kpath.klabels)
    plot!(p_mod, size=(800,500))
    ylims!(p_mod, define_ylims(band_data_mod, εF))
    ylabel!(p_mod, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f\;(\mathrm{hartree})")
    title!(p_mod,"Modified kinetic term")

    # Save plots if a directory is provided
    if !isempty(savedir)
        savefig(p_ref, joinpath(savedir, "band_plot_ref.pdf"))
        savefig(p_std, joinpath(savedir, "band_plot_std.pdf"))
        savefig(p_mod, joinpath(savedir, "band_plot_mod.pdf"))
    end

    p_ref, p_std, p_mod
end

function compute_axis_attributes(path_section, num_k, ref_data)
    # k point as 1D axis coordinates
    tmp_plot = DFTK.plot_band_data(ref_data.band_data; ref_data.scfres.εF,
                                   ref_data.kpath.klabels)
    ticks = only(xticks(tmp_plot))

    # Starting and end point coordinates
    k_start = ticks[1][findfirst(x->contains(x, path_section[1]), ticks[2])]
    k_end = ticks[1][findfirst(x->contains(x, path_section[2]), ticks[2])]
    # Axis
    x_axis = LinRange(k_start, k_end, num_k)

    k_start, k_end, x_axis    
end

"""
Plot band computed with the "focus_on_band" routine.
zoom_data contains values of the band for reference, standard and modified
kinetic term as well as the section of the path on which to plot.
"""
function plot_band(band_ref, band_std, bands_mod;
                   ref_data, savedir="",
                   i_derivative=zero(Int64))
    @assert (2+i_derivative ≤ length(band_ref.data)) "Order of derivative to high"

    # Extract bands and derivatives
    εn_ref = band_ref.data[2+i_derivative]
    εn_std = band_std.data[2+i_derivative]

    # Compute x_axis attributes
    path_section = band_std.path_section
    num_k = length(εn_ref)
    k_start, k_end, x_axis = compute_axis_attributes(path_section, num_k, ref_data)
    
    # Derivative symbol
    dev_symbl = ["", "∂", "∂^{2}"][i_derivative+1]

    # Plot band
    p = plot()
    (i_derivative == 0) && (plot!(p, x_axis, εn_ref,
                                  label=latexstring("\\varepsilon_{$(n)k}"),
                                  linecolor=:black, linestyle=:dash, linewidth=1.2,
                                  # markershape=:auto, markercolor=:match))
                                  )
                            )
    (i_derivative ≤ 1) && (plot!(p, x_axis, εn_std,
                                 label=latexstring(dev_symbl*"\\varepsilon_{$(n)k}^{E_c}"),
                                 linewidth=1.2)#, markershape=:auto, markercolor=:match)
                           )
    
    for (k, band_mod) in enumerate(bands_mod)
        εn_mod = band_mod.data[2+i_derivative]
        blow_up_rate = extract_blow_up_rate(band_mod.data[1][1])
        plot!(p, x_axis, εn_mod, label=latexstring(dev_symbl*
                     "\\tilde{\\varepsilon}_{$(n)k}^{E_c},\\; |⋅|^{$(blow_up_rate)}"),
              linewidth=1.2, linecolor=palette([:green, :red], length(bands_mod))[k],
              # markershape=:auto, markercolor=:match)
              )
    end
    
    plot!(p, size=(800,500), legend=:bottomleft)
    ylabel!(p, "eigenvalues (hartree)")
    xlabel!(p, "wave vector")
    xticks!(p, [k_start, k_end], path_section)
    !isempty(savedir) && (savefig(p,
                           joinpath(savedir, "band_$(n)_dev_$(i_derivative)_"*
                                    "$(path_section[1])_$(path_section[2]).pdf"))
                          )                          
    p
end
