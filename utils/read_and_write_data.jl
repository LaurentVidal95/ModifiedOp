function save_band_data(band_data, savedir::String)
    @assert !(isdir(savedir)) "$(savedir) already exists" 
    mkdir(savedir)
    @info "Saving band_data in directory $(savedir)"
    for ik in 1:length(band_data.λ)
        writedlm(joinpath(savedir,"X_$(ik).dat"), band_data.X[ik], ",")
        writedlm(joinpath(savedir,"λ_$(ik).dat"), band_data.λ[ik], ",")
    end
    nothing
end

function read_band_data(bs_basis, savedir::String)
    num_k = length(bs_basis.kpoints)
    X = [readdlm(joinpath(savedir,"X_$(ik).dat"), ',', ComplexF64) for ik in 1:num_k]
    λ = [readdlm(joinpath(savedir,"λ_$(ik).dat"), ',', Float64) for ik in 1:num_k]    
    εn = n->[εnk[n] for εnk in λ]
    ∂εn = n-> bands_derivative(εn(n), bs_basis.kcoords_global)
    (;bs_basis, X, λ), εn, ∂εn
end
function read_band_data(basis, kcoords, savedir::String)
    ksymop = [[one(SymOp)] for _ in 1:length(kcoords)]
    read_band_data(PlaneWaveBasis(basis, kcoords, ksymop), savedir)
end
