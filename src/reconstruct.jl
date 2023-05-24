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

function spectralfunction(chains::Vector{WilsonChain}; zaveraging=false, save=true, takezeroth=true)
    Nz = length(chains)
    freqs = Array{Vector{ComplexF64}}(undef, Nz)
    weights = Array{Vector{ComplexF64}}(undef, Nz)
    
    for (i, chain) in enumerate(chains)
        println("Diagonalising for (i, Nz) = ($(i), $(Nz))")
        freqs[i], weights[i] = spectralfunction(chain; takezeroth=takezeroth)
        if save == true
            mkpath("$(i)")
            open("$(i)/spec_f0.dat", "w") do f
                writedlm(f, hcat(real.(freqs[i]), real.(weights[i])))
            end
        end
    end

    if zaveraging == true
        freqs = (1/Nz) .* reduce(+, freqs)
        weights = (1/Nz) .* reduce(+, weights)
        if save == true
            open("spec_f0_zaveraged.dat", "w") do f
                writedlm(f, hcat(real.(freqs), real.(weights)))
            end
        end
    end
    
    return freqs, weights 
end

function spectralfunction(chain::WilsonChain; takezeroth=true)
    E = chain.E
    T = chain.T
    N = numberchannels(chain)
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

    if takezeroth == true
        H = Matrix(BlockTridiagonal(vecT, pushfirst!(vecE, zeros(N, N)), vecT))
        vals, vecs = eigen(H)
        weights = zeros(size(vecs, 2))
        for (i, k) in enumerate(eachcol(vecs))
            weights[i] = abs(k[1])^2
        end
        return vals, weights
    else
        H = Matrix(BlockTridiagonal(vecT[2:end], vecE, vecT[2:end]))
        vals, vecs = eigen(H)
        weights = zeros(size(vecs, 2))
        for (i, k) in enumerate(eachcol(vecs))
            weights[i] = abs(k[1])^2
        end
        return vals, weights
    end
end