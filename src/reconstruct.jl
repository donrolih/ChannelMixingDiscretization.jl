"""
    Given a star Hamiltonian reconstruct the original hybridization function.
"""
function reconstructhybri(star::StarHamiltonian, ωs; smear=0.1)
    Es = star.E
    Ts = star.T
    Elist = vcat(reverse(Es[-1]; dims=1), Es[1])
    Tlist = vcat(reverse(Ts[-1]; dims=1), Ts[1])
    if ndims(Elist) == 2
        # single band case
        println("Single band case")
        J, Nz = size(Elist)
        hybri = zeros(ComplexF64, length(ωs))
        println("Starting reconstruction of the hybridisation function!")
        println("Reconstructing ...")
        start = time()
        for (i, ω) in enumerate(ωs)
            Σ = complex(0.)
            for z in 1:Nz
                for j in 1:J
                    ti = Tlist[j, z]
                    ei = Elist[j, z]
                    gi = (ω + im*(smear/Nz)*max(1e-4, abs(ω))) - ei
                    Σ += conj(ti)*(1/gi)*ti
                end
            end
            hybri[i] = (im/(2*pi))*(Σ - Σ')
        end
        stop = time()
        elapsed = round(stop - start; digits=2)
        println("Done in $(elapsed) seconds!")
        return hybri/Nz
    else
        # multiband case
        println("Multi band case")
        J, Nz, Nband, _ = size(Elist)
        hybri = zeros(ComplexF64, length(ωs), Nband, Nband)
        println("Starting reconstruction of the hybridisation function!")
        println("Reconstructing ...")
        start = time()
        for (i, ω) in enumerate(ωs)
            Σ = zeros(ComplexF64, Nband, Nband)
            for z in 1:Nz
                for j in 1:J
                    ti = Tlist[j, z, :, :]
                    ei = Elist[j, z, :, :]
                    gi = (ω + im*(smear/Nz)*max(1e-4, abs(ω)))*I - ei
                    invgi = inv(gi)
                    Σ += ti'*invgi*ti
                end
            end
            hybri[i, :, :] = (im/(2*pi))*(Σ - Σ')
        end
        stop = time()
        elapsed = round(stop - start; digits=2)
        println("Done in $(elapsed) seconds!")
        return hybri/Nz
    end
end

# this reconstruction method is described in Appendix 2 of PRB 93, 035102 (2016)
"""
    Given a vector of chain Hamiltonians reconstruct the original hybridization function.
"""
function reconstructhybri(chains::Vector{WilsonChain}, ωs; smear=0.1)
    Nz = length(chains)
    # Getting Nbands is a bit clumsy, maybe I should write a function for this?
    Nbands = size(chains[1].E, 2)
    hybri = zeros(ComplexF64, length(ωs), Nbands, Nbands)
    start = time()
    println("Starting reconstruction of the spectrum with Nz = $(Nz)!")
    for (i, chain) in enumerate(chains)
        z = chain.z
        println("($(i)/$(Nz)): reconstructing spectrum for z = $(z) ...")
        E = chain.E
        T = chain.T
        println(size(E))
        println(size(T))
        for (j, ω) in enumerate(ωs)
            Σ = zeros(ComplexF64, Nbands, Nbands)
            # we start at the last site of the chain (reverse is necessary!)
            for k in size(E, 1):-1:1
                e = E[k, :, :]
                t = T[k, :, :]
                invG = (ω + im*((abs(ω) + 1e-10)*(smear/Nz)))*I - e - Σ
                G = inv(invG)
                Σ = t'*G*t
            end
            hybri[j, :, :] += (im/(2*pi))*(Σ - Σ')
        end
    end
    stop = time()
    elapsed = round(stop - start; digits=2)
    println("Done in $(elapsed) seconds!")
    # average over all twisting parameters
    return hybri/Nz
end

function spectralfunction(chains::Vector{WilsonChain}; save=true)
    Nz = length(chains)
    freqs = Array{Vector{ComplexF64}}(undef, Nz)
    weights = Array{Vector{ComplexF64}}(undef, Nz)
    
    for (i, chain) in enumerate(chains)
        println("Diagonalising for (i, Nz) = ($(i), $(Nz))")
        freqs[i], weights[i] = spectralfunction(chain)
        if save == true
            mkpath("spectral_functions")
            open("spectral_functions/$(i)_spec_f0.dat", "w") do f
                writedlm(f, hcat(real.(freqs[i]), real.(weights[i])))
            end
        end
    end
    
    return freqs, weights 
end

function buildhamiltonian(chain::WilsonChain; with_imp=false)
    E = chain.E
    T = chain.T

    J = numbersites(chain)
    vecE = Array{Array{ComplexF64, 2}, 1}(undef, J)
    vecT = Array{Array{ComplexF64, 2}, 1}(undef, J)
    vecT_adj = Array{Array{ComplexF64, 2}, 1}(undef, J)
    for j in axes(T, 1)
        t = T[j, :, :]
        vecT[j] = t
        vecT_adj[j] = t'
        vecE[j] = E[j, :, :]
    end
    
    if with_imp
        # with the impurity site
        H = Matrix(BlockTridiagonal(vecT[1:end], vcat([zeros(ComplexF64, 2, 2)], vecE), vecT_adj[1:end]))
        return H
    else
        # remove the impurity site
        H = Matrix(BlockTridiagonal(vecT[2:end], vecE, vecT_adj[2:end]))
        return H
    end
end

function buildhamiltonian(star::StarHamiltonian; znum=1)
    Epos = convert.(ComplexF64, star.E[1][:, znum, :, :])
    Eneg = convert.(ComplexF64, star.E[-1][:, znum, :, :])
    T = convert.(ComplexF64, getTlist(star)[:, znum, :, :])
    J, N, _ = size(Epos)
    H = zeros(ComplexF64, (2J + 1)*N, (2J + 1)*N)
    vecEpos = Array{Array{ComplexF64, 2}, 1}(undef, J)
    vecEneg = Array{Array{ComplexF64, 2}, 1}(undef, J)
    for j in 1:J 
        vecEpos[j] = Epos[j, :, :]
        vecEneg[j] = Eneg[j, :, :]
    end
    vecEneg = reverse(vecEneg)
    H = Matrix(BlockDiagonal(vcat([zeros(ComplexF64, N, N)], vecEneg, vecEpos)))
    v1 = reshape(transpose(T[:, 1, :]), 2J*N)
    v2 = reshape(transpose(T[:, 2, :]), 2J*N)
    H[1, (N + 1):end] = v1
    H[2, (N + 1):end] = v2
    H[(N + 1):end, 1] = conj.(v1)
    H[(N + 1):end, 2] = conj.(v2)
    return H
end

function spectralfunction(chain::WilsonChain)
    H = buildhamiltonian(chain)
    vals, vecs = eigen(H)
    weights = zeros(size(vecs, 2))
    for (i, k) in enumerate(eachcol(vecs))
        weights[i] = abs(k[1])^2
        # weights[i] = real.(k[1]'*k[2])
    end
    return vals, weights
end

################
# BROADENING   #
################

function gaussiankernel(ω, E, weight; η=0.05)
    σ = η*abs(E)
    A = (1/sqrt(2pi))*(1/σ)
    value = A*weight*exp(-(ω - E)^2/(2*σ^2))
    return value
end

function broaden(ωs, energies, weights)
    spectralfunction = zeros(length(ωs))
    for (i, E) in enumerate(energies)
        weight = weights[i]
        spectralfunction += gaussiankernel.(ωs, E, weight)
    end
    return spectralfunction
end

function broadenaverage(lb, ub, Nω, Nz; save=true)
    ωs = range(lb, ub, Nω)
    result = zeros(Nω)
    
    for i in 1:Nz
        data = readdlm("spectral_functions/$(i)_spec_f0.dat")
        energies = data[:, 1]
        weights = data[:, 2]
        result += broaden(ωs, energies, weights)
    end
    result = result ./ Nz

    open("spectral_function.dat", "w") do f
        writedlm(f, [ωs result])
    end
    return ωs, result
end