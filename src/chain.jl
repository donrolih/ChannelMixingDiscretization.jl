struct WilsonChain
    E::Array{Complex{BigFloat}, 3}
    T::Array{Complex{BigFloat}, 3}
    z::Tuple{Int64, Int64}
end

function numberchannels(chain::WilsonChain)
    return size(chain.E, 2)
end

function numbersites(chain::WilsonChain)
    return size(chain.E, 1)
end

function getTlist(star::StarHamiltonian)
    return vcat(reverse(star.T[-1]; dims=1), star.T[1])
end

function getElist(star::StarHamiltonian)
    return vcat(reverse(star.E[-1]; dims=1), star.E[1])
end

function maptochains(starH::StarHamiltonian; m=nothing)
    Elist = getElist(starH)
    Tlist = getTlist(starH)
    zs = starH.zs
    # I probably have to again divide separate cases for Nbands = 1 and other
    # if Elist and Tlist are matrices (not tensors) then we have one band
    if isa(Elist, Matrix)
        J, Nz = size(Elist)
        Nbands = 1
        chains = Vector{WilsonChain}(undef, Nz)
        for (i, z) in enumerate(zs)
            println("Mapping to a Wilson chain for twisting parameter z = $(z) ...")
            start = time()
            ti = Tlist[:, i]
            # make a matrix
            ti = ti[:, :]
            F = qr(ti)
            H = diagm(Elist[:, i])
            q = Matrix(F.Q)
            diags = tridiagonalize(H, q, m)
            # construct E and T matrices along the chain
            E = diags[0]
            T = zeros(Complex{BigFloat}, 1, Nbands, Nbands)
            T[1, :, :] = F.R
            # sub-diagonal
            T = vcat(T, diags[-1])
            chains[i] = WilsonChain(E, T, (i, Nz))
            stop = time()
            elapsed = stop - start
            println("Mapping duration z = $(z): $(elapsed)")
        end
    else
        # we have a multiband case
        J, Nz, Nbands, _ = size(Elist)
        # store the data in a vector of dicts
        chains = Vector{WilsonChain}(undef, Nz)
        for (i, z) in enumerate(zs)
            println("Mapping to a Wilson chain for twisting parameter z = $(z) ...")
            start = time()
            ti = Tlist[:, i, :, :]
            ti = permutedims(reshape(permutedims(ti, [3, 2, 1]), (Nbands, :)), Nbands:-1:1)
            F = qr(ti)
            # Matrix(F.Q) returns the 'thin' Q-matrix with ortonormal columns
            # permuting this is necessary because Julia is a column-major language
            q = permutedims(reshape(permutedims(Matrix(F.Q), [2, 1]), (Nbands, :)), Nbands:-1:1)
            H = BlockDiagonal([Elist[k, i, :, :] for k in 1:J])
            # println(typeof(H))
            diags = tridiagonalize(H, q, m)
            # construct E and T matrices along the chain
            # the diagonal
            E = diags[0]
            T = zeros(Complex{BigFloat}, 1, Nbands, Nbands)
            T[1, :, :] = F.R
            # sub-diagonal
            T = vcat(T, diags[-1])
            chains[i] = WilsonChain(E, T, (i, Nz))
            stop = time()
            elapsed = stop - start
            println("Mapping duration z = $(z): $(elapsed)")
        end
    end
    return chains
end