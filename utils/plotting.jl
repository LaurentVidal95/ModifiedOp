# Default plot parameters
pyplot()
myblue=RGB(95/255,133/255,255/255)
mygreen=RGB(60/255,195/255,26/255)
default(fontfamily="serif",
        linewidth=0.75, framestyle=:box,
        label=nothing, grid=:true,
        linecolor=myblue,
        tickfontsize=14,
        titlefontsize=14,
        legendfontsize=14,
        guidefontsize=14)

"""
Plot gm blow up function on [0,1) given third part ha.
"""
function plot_blowup_function(blowup; plot_dir="")    
    x_axis = LinRange(0, 0.99, 10000)
    blowup_rate = blowup.p

    # Plot
    p = plot(x_axis, blowup.(x_axis), label=:none, xlims=[0,1])
    
    vline!(blowup.blowup_function.interval, label="smooth interpolation",
           linestyle=:dash, linecolor=mygreen)
    vline!([1], label=:none, linecolor=:black)
    plot!(x->x^2, label=L"x\mapsto x^2", linestyle=:dash, linecolor=:black)

    # Plot parameters
    plot!(legendfontsize=12, legend=:topleft)
    plot!(size=(750,500))
    xlabel!(L"x")
    ylabel!(L"\mathscr{G}(x)")
    title!("Blow up rate: "*latexstring("|\\cdot|^{-$(blowup_rate)}"))
    ylims!(0,5)

    # Save plot if needed
    !isempty(plot_dir) && (savefig(p, joinpath(plot_dir,"blowup_function.pdf")))

    p
end

"""
Plot M_EC along reference k-path
"""
function plot_M_EC(ref_data; plot_dir="")
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

    title!(split(ref_data.system.name,"_",keepempty=false)[1])
    xlabel!(L"\mathbf{k}"*"-point")
    ylabel!(L"M_{\mathrm{E}_\mathrm{c}}(\mathbf{k})-\mathrm{mean}(M_{\mathrm{E}_\mathrm{c}})")
    plot!(size=(500,500))
    xticks!(ticks...)

    !isempty(plot_dir) && (savefig(p, joinpath(plot_dir, "M_EC.pdf")))

    p
end

"""
Plot bandstructures with respective standard and modified kinetic terms
"""
function plot_bandstructures(plot_data; ref_data, plot_dir="")
    # Extract data
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

    # Needed by DFTK.
    occupation_threshold = DFTK.default_occupation_threshold()
    
    # Reference
    εF = DFTK.compute_occupation(band_data_ref.basis, band_data_ref.λ;
                                 occupation_threshold).εF
    λ_ref = [λk .- εF for λk in band_data_ref.λ]
    p_ref = DFTK.plot_band_data(band_data_ref; εF=εF, kpath.klabels,
                                linewitdh=1.2, linecolor=:black)
    ylims!(p_ref, define_ylims(band_data_ref, εF))
    plot!(p_ref, size=(500,500))
    ylabel!(p_ref, L"\varepsilon_{n,k}-\varepsilon_f\;(\mathrm{hartree})")
    xlabel!(p_ref," ")

    # Standard
    εF = DFTK.compute_occupation(band_data_std.basis, band_data_std.λ;
                                 occupation_threshold).εF
    λ_std = [λk .- εF for λk in band_data_std.λ]
    p_std = DFTK.plot_band_data(band_data_std; εF=εF, kpath.klabels)
    plot!(p_std, size=(500,500))
    ylims!(p_std, define_ylims(band_data_std, εF))
    xlabel!(p_std, " ")
    ylabel!(p_std, L"\varepsilon_{n,k}^{E_c}-\varepsilon_f\;(\mathrm{hartree})")
    title!(p_std,L"\mathrm{Standard\; kinetic\; term}")

    # Modified
    εF = DFTK.compute_occupation(band_data_mod.basis, band_data_mod.λ;
                                 occupation_threshold).εF
    λ_mod = [λk .- εF for λk in band_data_mod.λ]
    p_mod = DFTK.plot_band_data(band_data_mod; εF, kpath.klabels)
    plot!(p_mod, size=(500,500))
    ylims!(p_mod, define_ylims(band_data_mod, εF))
    ylabel!(p_mod, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f\;(\mathrm{hartree})")
    title!(p_mod,"Modified kinetic term")

    # Error plot
    num_bands = length(λ_ref[1])
    err_std = norm.(λ_ref - λ_std, 1) ./ num_bands
    err_mod = norm.(λ_ref - λ_mod, 1) ./ num_bands
    ticks=only(xticks(p_mod))
    x_axis = LinRange(ticks[1][1], ticks[1][end], length(λ_ref))
    p_err = plot(x_axis, err_std,
                 label=L"\varepsilon_{n,k}^{E_c}",
                 linecolor=myblue)
    plot!(p_err, x_axis, err_mod,
          label=L"\tilde{\varepsilon}_{n,k}^{E_c}",
          linecolor=mygreen)
    plot!(p_err, size=(500,500))
    xticks!(ticks...)
    xlabel!(L"\mathrm{wave\; vector\;}k ")
    ylabel!("Mean absolute deviation (hartree)")
          
    # Save plots if a directory is provided
    if !isempty(plot_dir)
        savefig(p_ref, joinpath(plot_dir, "band_plot_ref.pdf"))
        savefig(p_std, joinpath(plot_dir, "band_plot_std.pdf"))
        savefig(p_mod, joinpath(plot_dir, "band_plot_mod.pdf"))
        savefig(p_err, joinpath(plot_dir, "band_plot_errors.pdf"))
    end

    p_ref, p_std, p_mod, p_err
end

function merge_band_plots(p_ref, p_std, p_mod)
    # Extract plots data
    x_axis = [s[:x] for s in p_ref.series_list]
    y_ref = [s[:y] for s in p_ref.series_list]
    y_std = [s[:y] for s in p_std.series_list]
    y_mod = [s[:y] for s in p_mod.series_list]

    # Plot
    p = plot()
    plot!(size=(500, 500))
    xlabel!(p, L"\mathrm{wave\; vector\;}\mathbf{k}")
    ylabel!(p, "Fermi level-adjusted eigenvalues (Hartree)")

    # Split in two to avoid legend printing issues
    plot!(x_axis[1:end-1], y_ref[1:end-1], label=:none, linecolor=:black,
          linestyle=:dash, linewidth=1.2)
    plot!(x_axis[1:end-1], y_std[1:end-1], label=:none, linecolor=myblue)
    plot!(x_axis[1:end-1], y_mod[1:end-1], label=:none, linecolor=mygreen)

    plot!(x_axis[end], y_ref[end], label=L"\varepsilon_{n}(\mathbf{k})"*
          L"-\varepsilon\mathrm{F}(\varepsilon_{n})", linecolor=:black,
          linestyle=:dash, linewidth=1.2)
    plot!(x_axis[end], y_std[end], label=L"\varepsilon_{n}^{E_c}(\mathbf{k})"*
          L"-\varepsilon\mathrm{F}(\varepsilon_n^{E_c})", linecolor=myblue)
    plot!(x_axis[end], y_mod[end], label=L"\tilde{\varepsilon}_{n}^{E_c}(\mathbf{k})"*
          L"-\varepsilon\mathrm{F}(\tilde{\varepsilon}_{n}^{E_c})", linecolor=mygreen)

    xticks!(p, only(xticks(p_ref))...)
    p
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
                   ref_data, plot_dir="",
                   i_derivative=zero(Int64),
                   )
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

    # Plots
    p = plot()
    default(tickfontsize=17, titlefontsize=17, legendfontsize=17, guidefontsize=17)

    # Reference band
    plot!(p, x_axis, εn_ref,
          label=latexstring(dev_symbl*"\\varepsilon_{$(n)\\mathbf{k}}"),
          linecolor=:black, linestyle=:dash, linewidth=1.2,
          )

    # Standard band
    (i_derivative == 0) && (plot!(p, x_axis, εn_std,
                                 label=latexstring(dev_symbl*"\\varepsilon_{$(n)\\mathbf{k}}^{\\mathrm{E}_\\mathrm{c}}"),
                                 linewidth=1.2)
                           )
    # Modified band
    subset(tab, n) = [x for (i,x) in enumerate(tab) if rem(i,n)==0]
    for (k, band_mod) in enumerate(bands_mod[1:end])
        εn_mod = band_mod.data[2+i_derivative]
        blowup_rate = ["1/2", "3/2", "5/2"][k] #extract_blowup_rate(band_mod.data[1][1])

        # Avoid numerical instabilities due to FD derivative computation
        debug = div(length(εn_ref), length(εn_mod))
        x_axis_tmp = subset(x_axis, debug)

        plot!(p, x_axis_tmp, εn_mod, label=latexstring(dev_symbl*
                     "\\tilde{\\varepsilon}_{$(n)\\mathbf{k}}^{\\mathrm{E}_\\mathrm{c}}"*
                     ",\\; p=$(blowup_rate)}"),
              linewidth=1.2, linecolor=palette([:green, :red], length(bands_mod))[k],
              )
    end

    # Plot args
    plot!(p, size=(500,500), legend=:bottomleft)
    ylabel!(p, "Eigenvalues (hartree)")
    xlabel!(p, L"\mathbf{k}"*"-point")
    xticks!(p, [k_start, k_end], path_section)
    !isempty(plot_dir) && (savefig(p,
                           joinpath(plot_dir, "band_$(n)_dev_$(i_derivative)_"*
                                    "$(path_section[1])_$(path_section[2]).pdf"))
                          )                          
    p, (x_axis, k_start, k_end)
end
