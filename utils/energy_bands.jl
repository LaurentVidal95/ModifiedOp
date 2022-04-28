using StaticArrays

"""
Construct a PlaneWaveBasis given initial basis and
custom k-points cartesian coordinates.
"""
basis_given_kcoords(basis, kcoords) = PlaneWaveBasis(basis, kcoords,
                               [[one(SymOp)] for _ in 1:length(kcoords)])

"""
Generate k-points cartesian coordinates given starting point, end point,
and the number of wanted k-points.
"""
function generate_kpath(k_start, k_end, num_k)
    [SVector{3,Float64}(kcoord) for kcoord in
               map(x-> (1-x)*k_start .+ x*k_end, LinRange(0, 1, num_k))]    
end

"""
Compute finite difference derivative of the band εn.
εn is the vector of eigenvalues for each k-points.
"""
function band_derivative(εn::Vector{T}, kcoords) where {T<:Real}
    δk = norm(kcoords[2] .- kcoords[1])
    [(εn[i+1]-εn[i])/δk for i in 1:length(εn)-1]
end
band_derivative(εn, k_start, k_end, num_k) =
    bands_derivative(εn,  discretize_kpath(k_start, k_end, num_k))

"""
Compute the largest variation of ∂εn and returns a kpath with "num_k_out" points
around said variation in order to numercically zoom on it.
"""
function kpath_near_band_irregularity(∂εn, kcoords, num_k_out)
    # Compute point of larger derivative variation
    num_k_in = length(∂εn)
    jumps = [∂εn[i+1] - ∂εn[i] for i in 1:num_k_in-1]
    k_irr = findmax(jumps)[2]
    # Handle extreme points
    width_left = (k_irr<50) ? zero(Int64) : 50
    width_right = (k_irr>num_k_in-50) ? zero(Int64) : 50
    # Define starting and end point of path
    k_start = kcoords[k_irr - width_left]
    k_end = kcoords[k_irr + width_right]
    # Return kpath
    generate_kpath(k_start, k_end, num_k_out)
end
