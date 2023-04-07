struct WilsonChain{S<:Real}
    E::Array{Complex{S}, 3}
    T::Array{Complex{S}, 3}
    z::Float64
end

function maptochains(starH::StarModel; m=nothing)
    Elist = getElist(starH)
    Tlist = getTlist(starH)
    zs = starH.zs
    # I probably have to again divide separate cases for Nbands = 1 and other
    J, Nz, Nbands, _ = size(Elist)
    # store the data in a vector of dicts
    chains = Vector{WilsonChain}(undef, Nz)
    for (i, z) in enumerate(zs)
        ti = Tlist[:, i, :, :]
        ti = permutedims(reshape(permutedims(ti, [3, 2, 1]), (Nbands, :)), Nbands:-1:1)
        F = qr(ti)
        # Matrix(F.Q) returns the 'thin' Q-matrix with ortonormal columns
        q = permutedims(reshape(permutedims(Matrix(F.Q), [2, 1]), (Nbands, :)), Nbands:-1:1)
        println("Mapping to a Wilson chain for twisting parameter z = $(z) ...")
        H = BlockDiagonal([Elist[k, i, :, :] for k in 1:J])
        # println(typeof(H))
        diags = tridiagonalize(H, q, m)
        # construct E and T matrices along the chain
        # the diagonal
        E = diags[0]
        T = zeros(ComplexF64, 1, Nbands, Nbands)
        T[1, :, :] = F.R
        # sub-diagonal
        T = vcat(T, diags[-1])
        chains[i] = WilsonChain(E, T, z)
    end

    return chains
end