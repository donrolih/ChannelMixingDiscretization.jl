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
                # invG = (ω + im*(smear/Nz)*max(1e-4, abs(ω)))*I - e - Σ
                # invG = (ω + im*(smear/Nz)*maximum(abs, t))*I - e - Σ
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