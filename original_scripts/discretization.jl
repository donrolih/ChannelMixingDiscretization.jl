include("integrators.jl")
using DataInterpolations
using ChannelMixingDiscretization
using BlockDiagonals
using LinearAlgebra

struct AdaptiveGrid
    interp_iR::DataInterpolations.AbstractInterpolation
    Rtot::Real
    function AdaptiveGrid(ωs, weights)
        R = cumintegrate(ωs, weights)
        Rtot = R[end]
        interp_iR = LinearInterpolation(ωs, R)
        return new(interp_iR, Rtot)
    end
end

struct FixedGrid
    Λ::Real
    D::Real
    Δ::Real
end

FixedGrid() = FixedGrid(big"2.", big"1.", big"0.0")

function makeguiding(grid::AdaptiveGrid)
    iR = grid.interp_iR
    Rtot = grid.Rtot
    # WARNING: This only works for bandwidth 1!!!
    # 2024-02-16: Now fixed using iR(Rtot) = Dmax
    # Should I use a field in the struct AdaptiveGrid to store Dmax?
    return x -> x ≤ big"2." ? iR(Rtot) : iR(Rtot*(big"2." ^ (big"2." - x)))
end

function makeguiding(grid::FixedGrid)
    return x -> x ≤ big"2." ? grid.D : (grid.D - grid.Δ)*grid.Λ^(big"2." - x) + grid.Δ
end

struct Discretizer
    ν
    interp_R
    interp_iR
    interp_Rint
    bsign
    function Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype, xmax, Nx)
        grid = gridtype == :fixed ? FixedGrid() : AdaptiveGrid(bsign .* ωbranch, weightbranch)
        ν = makeguiding(grid)
        R = cumintegrate(bsign .* ωbranch, ρbranch)
        interp_R = LinearInterpolation(R, bsign .* ωbranch)
        interp_iR = LinearInterpolation(bsign .* ωbranch, R)
        if typeof(xmax) == BigFloat
            xs = range(big"1.", xmax, Nx)
            Rint = cumintegrate(xs, interp_R(ν.(xs)))
            interp_Rint = LinearInterpolation(Rint, xs)
            return new(ν, interp_R, interp_iR, interp_Rint, bsign)
        else
            xs = range(1., xmax, Nx)
            Rint = cumintegrate(xs, interp_R(ν.(xs)))
            interp_Rint = LinearInterpolation(Rint, xs)
            return new(ν, interp_R, interp_iR, interp_Rint, bsign)
        end
    end
end

Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype) = Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype, big"200.", 10000)

function w(disc::Discretizer, x::AbstractFloat)
    ν = disc.ν
    interp_R = disc.interp_R
    upperbound = ν(x)
    lowerbound = ν(x + 1)
    return interp_R(upperbound) .- interp_R(lowerbound)
end

w(disc::Discretizer, xs) = [w(disc, x) for x in xs]

function ε(disc::Discretizer, x)
    interp_iR = disc.interp_iR
    interp_Rint = disc.interp_Rint
    bsign = disc.bsign
    return bsign .* interp_iR(interp_Rint(x .+ 1) .- interp_Rint(x))
end

ε(disc::Discretizer, xs::Vector{AbstractFloat}) = [ε(disc, x) for x in xs]

function calculatediscretizers(ωs, ρeval, weights, bsign, gridtype)
    @assert (bsign == 1) | (bsign == -1) "the value of sign is not 1 or -1"
    mask = bsign == 1 ? ωs .≥ big"0." : ωs .< big"0."
    ωbranch = bsign == 1 ? ωs[mask] : reverse(ωs[mask])
    weightbranch = bsign == 1 ? weights[mask] : reverse(weights[mask])
    # ρeval is a vector or a matrix
    if isa(ρeval, Vector)
        # one band case
        rhovec = ρeval
        ρbranch = bsign == 1 ? rhovec[mask] : reverse(rhovec[mask])
        return [Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype)]
    else
        # multiband case
        _, Nbands = size(ρeval)
        discs = Vector{Discretizer}(undef, Nbands)
        for i in 1:Nbands
            rhovec = ρeval[:, i]
            ρbranch = bsign == 1 ? rhovec[mask] : reverse(rhovec[mask])
            discs[i] = Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype)
        end
        return discs
    end
end

function evaluatecoefficients(discpos::Vector{Discretizer}, discneg::Vector{Discretizer}, ρ, zs, J)
    Ts = Dict()
    Es = Dict()
    Nz = length(zs)
    Nbands = length(discpos)
    if Nbands == 1
        # one band case
        for bsign in [-1, 1]
            disc = bsign == 1 ? discpos : discneg
            E = zeros(Complex{BigFloat}, J, length(zs))
            T = zeros(Complex{BigFloat}, J, length(zs))
            for j in 1:J
                for (k, z) in enumerate(zs)
                    x = big(j + z)
                    ϵ = ε(disc[1], x)
                    # println(ϵ)
                    t = sqrt(w(disc[1], x))
                    E[j, k] = ϵ
                    T[j, k] = t
                end
            end
            Es[bsign] = E
            Ts[bsign] = T
        end
    else
        # multiband case
        for bsign in [-1, 1]
            disc = bsign == 1 ? discpos : discneg
            E = zeros(Complex{BigFloat}, J, Nz, Nbands, Nbands)
            T = zeros(Complex{BigFloat}, J, Nz, Nbands, Nbands)
            for j in 1:J
                for (k, z) in enumerate(zs)
                    # println("sign: ", bsign)
                    # println("j: ", j)
                    x = big(j + z)
                    ϵ = [ε(disc[i], x) for i in 1:Nbands]
                    ϵ = diagm(ϵ)
                    E[j, k, :, :] = ϵ
                    ϵiold = Inf
                    Ux = zeros(Complex{BigFloat}, Nbands)
                    for i in 1:Nbands
                        ϵi = ϵ[i, i]
                        if abs(ϵiold - ϵi) > 1e-100
                            M = Hermitian(ρ(ϵi))
                            # println("ϵi: ", ϵi)
                            # println("hybridization: ", M)
                            Ui = eigvecs(M)[:, i]
                            ϵiold = ϵi
                        else
                            println("Found a degeneracy at index J = $(J)!")
                            if i == 1 error("impossible!") end
                            # we have to choose the same vector as before
                            Ui = Ux[:, end]
                        end
                        Ux = hcat(Ux, Ui)
                    end
                    Ux = Ux[:, 2:end]
                    # println(Ux)
                    Td = zeros(Complex{BigFloat}, Nbands)
                    for i in 1:Nbands
                        Tdi = sqrt(w(disc[i], x)) .* Ux[:, i]
                        # println("ti", sqrt(w(disc[i], x)))
                        Td = hcat(Td, Tdi)
                    end
                    T[j, k, :, :] = Td[:, 2:end]'
                end
            end
            Es[bsign] = E
            Ts[bsign] = T
        end
    end
    return Ts, Es
end

function maptochainsSPSU2(starH::StarHamiltonian; m=nothing)
    Elist = ChannelMixingDiscretization.getElist(starH)
    Tlist = ChannelMixingDiscretization.getTlist(starH)
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
            diags = ChannelMixingDiscretization.tridiagonalize(H, q, m)
            # construct E and T matrices along the chain
            E = diags[0]
            T = zeros(Complex{BigFloat}, 1, Nbands, Nbands)
            T[1, :, :] = F.R
            # sub-diagonal
            T = vcat(T, diags[-1])
            chains[i] = WilsonChain(E, T, (i, Nz))
            stop = time()
            elapsed = stop - start
            println("Mapping duration for z = $(z): $(round(elapsed;  digits=2)) s")
        end
    else
        # we have a multiband case
        # here J is actually 2*number of j 
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
            # the sign has to fixed by hand (because it is arbitrary but we want the whole thing to be SU2 symmetric)!!
            # the sign is arbitrary because of QR decomposition! 
            # Equivalently, we could have used the square root of a positive definite function. There we would also have the liberty to choose either the positive or the negative branch.
            S = diagm([sign(real(F.R[1, 1])) < 0 ? -1. : 1, sign(real(F.R[2, 2])) < 0 ? 1. : -1.])
            R = S*F.R
            q = Matrix(F.Q)*S
            # println(q)
            H = BlockDiagonal([Elist[k, i, :, :] for k in 1:J])
            # println(typeof(H))
            diags = tridiagonalizeSPSU2(H, q, m)
            # construct E and T matrices along the chain
            # the diagonal
            E = diags[0]
            T = zeros(Complex{BigFloat}, 1, Nbands, Nbands)
            T[1, :, :] = R
            # sub-diagonal
            T = vcat(T, diags[-1])
            chains[i] = WilsonChain(E, T, (i, Nz))
            stop = time()
            elapsed = stop - start
            println("Mapping duration for z = $(z): $(round(elapsed;  digits=2)) s")
        end
    end
    return chains
end

function maptochainsONEBAND(starH::StarHamiltonian; m=nothing)
    Elist = ChannelMixingDiscretization.getElist(starH)
    Tlist = ChannelMixingDiscretization.getTlist(starH)
    zs = starH.zs
    # I probably have to again divide separate cases for Nbands = 1 and other
    # if Elist and Tlist are matrices (not tensors) then we have one band
    if isa(Elist, Matrix)
        J, Nz = size(Elist)
        Nbands = 1
        chains = Vector{WilsonChain}(undef, Nz)
        for (i, z) in enumerate(zs)
            @info "Mapping to a Wilson chain for twisting parameter z = $(z) ..."
            start = time()
            ti = Tlist[:, i]
            # make a matrix
            ti = ti[:, :]
            F = qr(ti)
            H = diagm(Elist[:, i])
            # we want F.R to be positive
            # thus we have to multiply it by its sign
            S = sign(real(F.R[1, 1]))
            R = S*F.R
            q = Matrix(F.Q)*S
            diags = tridiagonalizeONEBAND(H, q, m)
            # construct E and T matrices along the chain
            E = diags[0]
            T = zeros(Complex{BigFloat}, 1, Nbands, Nbands)
            T[1, :, :] = R
            # sub-diagonal
            T = vcat(T, diags[-1])
            chains[i] = WilsonChain(E, T, (i, Nz))
            stop = time()
            elapsed = stop - start
            @info "Mapping duration for z = $(z): $(round(elapsed;  digits=2)) s"
        end
    else
        error("We should have a one band case!")
    end
    return chains
end

function tridiagonalizeSPSU2(A, q, m)
    n = size(q, 2)
    # set the number of iterations as 2J
    m = m === nothing ? Int(size(A)[1]/n) : m
    # Ss = zeros(Float64, m, n, n)
    Al = zeros(Complex{BigFloat}, m, n, n)
    B = zeros(Complex{BigFloat}, m+1, n, n)
    Bu = zeros(Complex{BigFloat}, m, n, n)
    Q = hcat(zeros(BigFloat, size(q)), q)
    for j in 1:m
        Qj = Q[:, end-(n-1):end]
        Uj = A*Qj - Q[:, (end - 2n + 1):(end - n)]*(B[j, :, :]')
        Aj = Qj'*Uj
        Rj = Uj - Qj*Aj
        # ortho = dot(Rj[:, 1], Rj[:, 2])
        # println("are the columns of Rj orthogonal? ", norm(ortho) < 1e-20 ? true : false, " inner product: ", ortho)
        # println("is Rj*Rj' positive definite? Eigenvalues: ", eigvals(Rj'*Rj))
        # Answer: yes!!!
        # Then we don't need to do QR! We basically just have to normalize the columns of Rj, because off diagonal elements in the R matrix are proportional to the inner product of columns of Rj
        F = qr(Rj)
        # new Q
        # S = diagm(sign.(real.(diag(F.R))))
        # S = I
        # we have to choose S s.t. the symmetry properties are satisfied
        S = diagm([sign(real(F.R[1, 1])) < 0 ? -1. : 1., sign(real(F.R[2, 2])) < 0 ? 1. : -1.])
        # Ss[j, :, :] = S
        R = S*F.R
        newQj = Matrix(F.Q)*S
        # println(typeof(newQj))
        # # println(size(Rj))
        # Rjnorms = zeros(Complex{BigFloat}, size(Rj, 2))
        # for (i, col) in enumerate(eachcol(Rj))
        #     colnorm = norm(col)
        #     Rjnorms[i] = colnorm
        # end
        # R = diagm(Rjnorms)
        # foreach(normalize!, eachcol(Rj))
        # newQj = Rj
        # reorthogonalize
        newQj = ChannelMixingDiscretization.gramschmidt(newQj, Q[:, n+1:end])
        # println("new Q block vector orthogonalization: ", newQj'*newQj)
        Q = hcat(Q, newQj)
        
        if (j != m) && norm(R, 1) < 1e-20
            println("WARNING! Bad Krylov space!")
        end
        
        # store the results
        Al[j, :, :] = Aj
        B[j + 1, :, :] = R
        Bu[j, :, :] = R'
    end
    Bl = B[2:end-1, :, :]
    # return a dictionary of arrays; keys denote the (off)-diagonal index
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end

function tridiagonalizeONEBAND(A, q, m)
    n = size(q, 2)
    # set the number of iterations as 2J
    m = m === nothing ? Int(size(A)[1]/n) : m
    # Ss = zeros(Float64, m, n, n)
    Al = zeros(Complex{BigFloat}, m, n, n)
    B = zeros(Complex{BigFloat}, m+1, n, n)
    Bu = zeros(Complex{BigFloat}, m, n, n)
    Q = hcat(zeros(BigFloat, size(q)), q)
    for j in 1:m
        Qj = Q[:, end-(n-1):end]
        Uj = A*Qj - Q[:, (end - 2n + 1):(end - n)]*(B[j, :, :]')
        Aj = Qj'*Uj
        Rj = Uj - Qj*Aj
        # ortho = dot(Rj[:, 1], Rj[:, 2])
        # println("are the columns of Rj orthogonal? ", norm(ortho) < 1e-20 ? true : false, " inner product: ", ortho)
        # println("is Rj*Rj' positive definite? Eigenvalues: ", eigvals(Rj'*Rj))
        # Answer: yes!!!
        # Then we don't need to do QR! We basically just have to normalize the columns of Rj, because off diagonal elements in the R matrix are proportional to the inner product of columns of Rj
        F = qr(Rj)
        # new Q
        # S = diagm(sign.(real.(diag(F.R))))
        # S = I
        # we have to choose S s.t. the symmetry properties are satisfied
        S = sign(real(F.R[1, 1]))
        # Ss[j, :, :] = S
        R = S*F.R
        newQj = Matrix(F.Q)*S
        # println(typeof(newQj))
        # # println(size(Rj))
        # Rjnorms = zeros(Complex{BigFloat}, size(Rj, 2))
        # for (i, col) in enumerate(eachcol(Rj))
        #     colnorm = norm(col)
        #     Rjnorms[i] = colnorm
        # end
        # R = diagm(Rjnorms)
        # foreach(normalize!, eachcol(Rj))
        # newQj = Rj
        # reorthogonalize
        newQj = ChannelMixingDiscretization.gramschmidt(newQj, Q[:, n+1:end])
        # println("new Q block vector orthogonalization: ", newQj'*newQj)
        Q = hcat(Q, newQj)
        
        if (j != m) && norm(R, 1) < 1e-20
            @warn "WARNING! Bad Krylov space!"
        end
        
        # store the results
        Al[j, :, :] = Aj
        B[j + 1, :, :] = R
        Bu[j, :, :] = R'
    end
    Bl = B[2:end-1, :, :]
    # return a dictionary of arrays; keys denote the (off)-diagonal index
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end