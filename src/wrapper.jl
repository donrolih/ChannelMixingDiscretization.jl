# Wrapper functions for the user

function discretize(ωs::Vector{BigFloat},
                    ρs::Vector{BigFloat},
                    mesh_min, mesh_max, mesh_ratio,
                    J, zs, gridtype,
                    Λ, gap, D;
                    savechain=true,
                    nrg_generatefolders=true,
                    )
    # this is single-channel case
    # the input is a tabulated function (ωs, ρs)
    
    interpolated_ρ = DataInterpolations.LinearInterpolation(ρs, ωs)
    
    freqs = logmesh(mesh_min, mesh_max, mesh_ratio, big"1e-100")
    
    # single-channel case, DOS values and weights  are the same as the input ρs
    values = interpolated_ρ.(freqs)
    weights = values

    gridparams = Dict{String, AbstractFloat}()
    gridparams["Lambda"] = Λ
    gridparams["gap"] = gap
    gridparams["halfbandwidth"] = D

    # calculate the discretizer objects
    posdiscs = calculatediscretizers(freqs, values, weights, 1, gridtype, gridparams)
    negdiscs = calculatediscretizers(freqs, values, weights, -1, gridtype, gridparams)
    
    # evaluate the star representation coeffients
    Es, Ts = evaluatecoefficients(posdiscs[1], negdiscs[1], J, zs)
    
    # build the star Hamiltonian
    starH = StarHamiltonian(Ts, Es, zs)

    # Wilson chain
    chains = maptochainsONEBAND(starH; m=2J)
    
    # save chain coeffients in a separate folder
    if savechain savechains(chains) end

    # generate NRG folders
    if nrg_generatefolders nrgfilesONEBAND(chains) end
    return starH, chains
end

function discretize(ωs::Vector{T},
                    ρs::Vector{Matrix{T}},
                    mesh_min, mesh_max, mesh_ratio, mesh_accumulation,
                    J, zs, gridtype, Λ, D;
                    savechain=true,
                    nrg_generatefolders=true
                    ) where T <: AbstractFloat
    # this is multi-channel case
    # the input is a tabulated function (ωs, ρs)
    # here, each element of ρs is a matrix of size (N, N)
    # where N is the number of channels

    # we need to interpolate each matrix element
    # i.e., we need ω -> ρ_ij(ω) for i,j = 1, ..., N
    # DataInterpolations package does this for us
    # interpolated_ρ(ω) = [ρ11(ω) ρ12(ω); ρ21(ω) ρ22(ω)]
    
    interpolated_ρ = DataInterpolations.LinearInterpolation(ρs, ωs)
    
    freqs = logmesh(mesh_min, mesh_max, mesh_ratio, mesh_accumulation)
    
    weights, ρeval = weightsDOS(interpolated_ρ.(freqs))

    gridparams = Dict{String, AbstractFloat}()
    gridparams["Lambda"] = Λ
    gridparams["gap"] = mesh_accumulation
    gridparams["halfbandwidth"] = D

    # discretization of the hybridization function
    posdiscs = calculatediscretizers(freqs, ρeval, weights, 1, gridtype, gridparams)
    negdiscs = calculatediscretizers(freqs, ρeval, weights, -1, gridtype, gridparams)

    # evaluate at points x = j + z
    Es, Ts = evaluatecoefficients(posdiscs, negdiscs, J, zs, interpolated_ρ)

    # collect the coeffients into star representation struct
    star = StarHamiltonian(Ts, Es, zs)

    # perform block Lanczos tridiagonalisation for each z number
    chains = maptochainsSPSU2(starH; m=2J)

    return star, chains
end

# Calculate the weight function and DOS for each channel
function weightsDOS(ρs::Vector{Matrix{T}}) where T <: AbstractFloat
    # get weights as a norm
    # Q: is there a better "measure" for the weights
    # Norm is λ₁² + λ₂²; is there a better way to define the weight?
    weights = [norm(ρ) for ρ in ρs]

    # diagonalize and get the eigenvalues
    ρeval = [eigvals(M) for M in ρs]
    ρeval = reduce(vcat, transpose.(ρeval))

    return weights, ρeval
end