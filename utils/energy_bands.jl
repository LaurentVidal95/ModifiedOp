function bands_derivative(εn, kcoords)
    δk = norm(kcoords[2] .- kcoords[1])
    [(εn[i+1]-εn[i])/δk for i in 1:length(εn)-1]
end
bands_derivative(εn, k_start, k_end, num_k) =
    bands_derivative(εn,  discretize_kpath(k_start, k_end, num_k))

function bands_irregularity(εn, k_coords)
    num_k = length(k_coords)
    ∂εn = bands_derivative(εn, k_coords)
    jumps = [∂εn[i+1] - ∂εn[i] for i in 1:num_k-2]
    k_mid = findmax(jumps)[2]
    k_start_zoom = k_coords[k_mid - 5]; k_end_zoom = k_coords[k_mid + 5];
    k_mid, (k_start_zoom, k_end_zoom)
end
