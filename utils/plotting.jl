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
function plot_gm(ha; savedir="")
    x_axis = LinRange(0, 0.99, 100)
    p = plot(x_axis, abs.(gm.(x_axis, ha)), label=:none, xlims=[0,1])
    vline!([1/2,3/4], label=L"C^2"*" spline interpolation",
           linestyle=:dash, linecolor=mygreen)
    vline!([1], label=:none, linestyle=:dash)
    plot!(x->x^2, label=L"x\mapstox^2", linestyle=:dash, linecolor=:black)
    plot!(legendfontsize=12, legend=:topleft)
    plot!(size=(750,500))
    ylims!(0,10)
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
function plot_band(zoom_data; ref_data, savedir="")
    # Extract bands and derivatives
    εn_ref = zoom_data.ref_data[2]
    εn_std = zoom_data.std_data[2]
    εn_mod = zoom_data.mod_data[2]

    ∂εn_std = zoom_data.std_data[3]
    ∂εn_mod = zoom_data.mod_data[3]
    
    # Compute x_axis attributes
    path_section = zoom_data.path_section
    num_k = length(εn_ref)
    k_start, k_end, x_axis = compute_axis_attributes(path_section, num_k, ref_data)

    # Plot band
    p = plot(x_axis, εn_ref, label=latexstring("\\varepsilon_{$(n)k}"),
             linecolor=:black, linestyle=:dash, linewidth=1.2)
    plot!(p, x_axis, εn_std, label=latexstring("\\varepsilon_{$(n)k}^{E_c}"),
          linewidth=1.2)
    plot!(p, x_axis, εn_mod, label=latexstring("\\tilde{\\varepsilon}_{$(n)k}^{E_c}"),
          linewidth=1.2, linecolor=mygreen)
    plot!(p, size=(800,500), legend=:topright)
    ylabel!(p, "eigenvalues (hartree)")
    xlabel!(p, "wave vector")
    xticks!(p, [k_start, k_end], path_section)

    # Plot first derivative
    p_std = plot(x_axis[1:end-1], ∂εn_std, label=latexstring("\\partial_k\\varepsilon_{$(n)}^{E_c}"))
    p_mod = plot(x_axis[1:end-1], ∂εn_mod, label=latexstring("\\partial_k\\tilde{\\varepsilon}_{$(n)}^{E_c}"))
    plot!(p_std, size=(800,500)); plot!(p_mod, size=(800,500))
    xlabel!(p_mod, "wave vector")
    # ylabel!(p_std, "FD derivative")
    # ylabel!(p_mod, "FD derivative")
    xticks!(p_mod, [k_start, k_end], path_section)

    p_dev = plot(p_std, p_mod, layout=(2,1))

    
    if !isempty(savedir)
        savefig(p, joinpath(savedir, "band_$(n)_$(path_section[1])_"*
                            "$(path_section[2]).pdf"))
        savefig(p_dev, joinpath(savedir, "band_$(n)_derivative_"*
                            "$(path_section[1])_$(path_section[2]).pdf"))
    end
    p, p_dev
end

function plot_derivatives_wr_blow_up_rate(n, blow_up_rates; ref_data,
                                          num_k=100, path_section=ref_data.kpath.kpath[1][1:2],
                                          savedir="")
    # Compute bands for modified kinetic terms for all given blow_up rates
    datas = focus_on_band.(n, blow_up_rates; ref_data, num_k, path_section, only_mod=true)

    # Compute x_axis attributes
    k_start, k_end, x_axis = compute_axis_attributes(path_section, num_k, ref_data)

    # Plot
    p1 = plot()
    p2 = plot()
    for (iε, ε) in enumerate(blow_up_rates)
        plot!(p1,  x_axis, datas[iε].mod_data[2],
              label=latexstring("|\\cdot|^{$(ε)}"),
              linecolor=iε)
        plot!(p2,  x_axis[1:end-1], datas[iε].mod_data[3],
              label=latexstring("|\\cdot|^{$(ε)}"),
              linecolor=iε)
    end
    plot!(p1, size=(800,500)); plot!(p2, size=(800,500))
    xlabel!(p2, "wave vector")
    xticks!(p2, [k_start, k_end], path_section)
    p12 = plot(p1, p2, layout=(2,1))

    !(isempty(savedir)) && (savefig(p12, joinpath(savedir, "band_$(n)_different_blowup_rates.pdf")))
    p12
end
