"""
Contains all routines to create and analyze a kpath with given
start point, end point and number of k_points.
"""
basis_given_kcoords(basis, kcoords) = PlaneWaveBasis(basis, kcoords,
                               [[one(SymOp)] for _ in 1:length(kcoords)])

function discretize_kpath(k_start, k_end, resolution)
    [SVector{3,Float64}(kcoord)
     for kcoord in map(x-> (1-x)*k_start .+ x*k_end, LinRange(0, 1, resolution))]
end

function compute_Gvec_jumps_along_path(basis)
    Xk_lengths = [length(G_vectors(basis, kpt)) for kpt in basis.kpoints]
    mean = round(Int64, sum(Xk_lengths) / length(Xk_lengths))
    @show num_k = length(Xk_lengths)
    jumps = [abs(Xk_lengths[i+1] - Xk_lengths[i]) for i in 1:num_k-1]
    (Xk_lengths .- mean), mean, jumps, findmax(jumps)
end

### Plots
function plot_info_kpath(basis_in, k_start, k_end, num_k)
    # Compute reduced kcoords and corresponding basis
    kcoords = discretize_kpath(k_start, k_end, num_k)
    num_k = length(kcoords)
    basis = basis_given_kcoords(basis_in, kcoords)

    # Plot jumps along path
    dist_to_mean_G, mean_G_vec, jumps, kpt_of_interest =
        compute_Gvec_jumps_along_path(basis)
    p = plot(collect(1:num_k)/num_k, dist_to_mean_G, markershape=:circle, label=:none)
    plot!(size=(750,500))
    hline!([0], linestyle=:dash, color=:black, label="mean = $(mean_G_vec)")
    plot!(legend=:topleft)
    xlabel!("k_path", xlabelfonsize=12)
    ylabel!("Distance to mean number of G vectors", ylabelfonsize=12)
    title!("Variations of number of G vectors along k path:"*
           "\n from $(round.(k_start; digits=3)) to $(round.(k_end; digits=3))")
    plot!(legendfontsize=12)

    # Return infos
    p, (kcoords, jumps, kpt_of_interest)
end
