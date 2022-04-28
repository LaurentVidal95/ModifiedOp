# Default plot parameters
pyplot()
default(fontfamily="serif",
        linewidth=0.75, framestyle=:box,
        label=nothing, grid=:true,
        linecolor=RGB(95/255,133/255,255/255),
        titlefontsize=22,
        tickfontsize=12,
        legendfontsize=12,
        guidefontsize=18)


"""
Plot gm blow up function on [0,1) given third part ha.
"""
function plot_gm(ha)
    x_axis = LinRange(0, 0.99, 100)
    p = plot(x_axis, abs.(gm.(x_axis, ha)), label=:none, xlims=[0,1])
    vline!([1/2,3/4], label="C^2 spline interpolation", linestyle=:dash)
    vline!([1], label=:none, linestyle=:dash)
    plot!(x->x^2, label="x->x^2", linestyle=:dash, color=:black)
    plot!(legendfontsize=12, legend=:topleft)
    plot!(size=(750,500))
    ylims!(0,10)
    p
end

"""
Plot M_EC along reference k-path
"""
function plot_M_EC(ref_data)
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
             
    xlabel!("wave vector")
    ylabel!(L"M_{E_c}-\mathrm{mean}(M_{E_c})")
    plot!(size=(800,500))
    xticks!(ticks...)
    p
end

"""
Plot bandstructures with respective standard and modified kinetic terms
"""
function bandplots(plot_data; ref_data, savedir="")
    # Extract data
    ρ_ref = ref_data.scfres.ρ
    band_data_std, λ_ref, ∂λ_ref = plot_data.std_data
    band_data, λ, ∂λ = plot_data.reg_data
    
    kpath = ref_data.kpath
    kcoords = kpath.kcoords
    num_k = length(kcoords)

    # Small function that computes the fermi level and min and max eigenvalue
    # so that the shift due to the change of kinetic term does not appear.

    function fermi_and_ylims(data, kcoords)
        εF = DFTK.fermi_level(data[1].basis, data[1][2])
        ymin = findmin(data[1][2])[1][1]
        ymax = findmax(data[1][2])[1][end]
        εF, (ymin, ymax)
    end

    # Reference
    p_ref = DFTK.plot_band_data(ref_data.band_data;
                                ref_data.scfres.εF, kpath.klabels)
    plot!(p_ref, size=(800,500))
    ylabel!(p_ref, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f")
    title!(p_ref,"Reference band diagram")

    # Standard
    εF, (ymin, ymax) = fermi_and_ylims(plot_data.std_data, kcoords)
    p_std = DFTK.plot_band_data(band_data_std; εF, kpath.klabels)
    plot!(p_std, size=(800,500))
    ylims!(p_std, (ymin, ymax))
    xlabel!(p_std, " ")
    ylabel!(p_std, L"\varepsilon_{n,k}^{E_c}-\varepsilon_f")
    title!(p_std,L"\mathrm{Standard\; kinetic\; term}")

    # Regularized
    εF, (ymin, ymax) = fermi_and_ylims(plot_data.reg_data, kcoords)
    p_reg = DFTK.plot_band_data(band_data; εF, kpath.klabels)
    plot!(p_reg, size=(800,500))
    ylims!(p_reg, (ymin, ymax))
    ylabel!(p_reg, L"\tilde{\varepsilon}_{n,k}^{E_c}-\varepsilon_f")
    title!(p_reg,"Modified kinetic term")

    # Save plots if a directory is provided
    if !isempty(savedir)
        savefig(p_ref, joinpath(savedir, "band_plot_ref.pdf"))
        savefig(p_std, joinpath(savedir, "band_plot_std.pdf"))
        savefig(p_reg, joinpath(savedir, "band_plot_reg.pdf"))
    end

    p_ref, p_std, p_reg
end

"""
Zoom on the largest irregularity. Plot ref, std and reg on the same plot.
BEWARE: for now mainly focus on crossings. Must find something to avoid
focusing on crossing. Add a user defined window of searching ?
"""
function zoom_on_band_irregularity(i_band, plot_data; ref_data, num_k=500)
    # Extract data from plot_data and ref_data
    ρ_ref = ref_data.scfres.ρ
    εn_std = plot_data.std_data[2]
    ∂εn_std = plot_data.std_data[3]
    basis_std = plot_data.std_data[1].basis
    basis_reg = plot_data.reg_data[1].basis
    
    # Define k-path in the neighbourhood of largest band irregularity with "num_k" points
    kcoords_zoom = kpath_near_band_irregularity(∂εn_std, num_k)
    
    # Compute bands around irregularity
    @info "Compute given band in the neighbourhood of"*
        " its biggest irregularity with $(num_k) points"
    # ref
    tmp_band_data = compute_bands(ref_data.scfres.basis,
                                  kcoords_zoom, n_bands=i_band+2, ρ=ρ_ref)
    εn_zoom_ref = [εnk[i_band] for εnk in tmp_band_data.λ]
    # std
    tmp_band_data = compute_bands(basis_std, kcoords_zoom, n_bands=i_band+2, ρ=ρ_ref)
    εn_zoom_std = [εnk[i_band] for εnk in tmp_band_data.λ]
    # reg
    tmp_band_data = compute_bands(basis_reg, kcoords_zoom, n_bands=i_band+2, ρ=ρ_ref)
    εn_zoom_reg = [εnk[i_band] for εnk in tmp_band_data.λ]

    # Plot
    p = plot(1:num_k, εn_zoom_ref, label=latexstring("\\varepsilon_{$(i_band)k}"),
             linecolor=:black, linestyle=:dash, linewidth=1.5)
    plot!(p, 1:num_k, εn_zoom_std, label=latexstring("\\varepsilon_{$(i_band)k}^{E_c}"),
          linewidth=1.5)
    plot!(p, 1:num_k, εn_zoom_reg, label=latexstring("\\tilde{\\varepsilon}_{$(i_band)k}^{E_c}"),
          linewidth=1.5, linecolor=RGB(60/255,195/255,26/255))
    plot!(p, size=(800,500), legend=:topright)
    ylabel!(p, "")
    xlabel!(p, "wave vector")
    xticks!(p, :none)
    p
end

"""
Compute the derivative of band number "n".
"""
function plot_band_derivative(n, plot_data; ref_data)
    # Compute ticks
    kpath = ref_data.kpath
    tmp_plot = DFTK.plot_band_data(ref_data.band_data; ref_data.scfres.εF, kpath.klabels)
    ticks = only(xticks(tmp_plot))

    # Extract band derivatives
    ∂εn_ref = band_derivative([εnk[n] for εnk in ref_data.band_data.λ], kpath.kcoords)
    ∂εn_std = plot_data.std_data[3]
    ∂εn_reg = plot_data.reg_data[3]

    # plot
    x_axis = LinRange(ticks[1][1], ticks[1][end], length(∂εn_ref))
    p_std = plot(x_axis, ∂εn_std[n], label=L"\partial_k\varepsilon_{n}^{E_c}")
    p_reg = plot(x_axis, ∂εn_reg[n], label=L"\partial_k\tilde{\varepsilon}_{n}^{E_c}")
    plot!(p_std, size=(800,500)); plot!(p_reg, size=(800,500))
    xlabel!(p_reg, "wave vector")
    xticks!(p_std, ticks...); xticks!(p_reg, ticks...);

    plot(p_std, p_reg, layout=(2,1))
end
