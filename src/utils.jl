########################
#GENERATING DATA
########################

"""
    Generate a vector of Nz twisting parameters, linearly distributed on (0, 1).
"""
function twistingparameters(Nz::Integer)
    return collect(range(0.5/Nz, 1 - 0.5/Nz, Nz))
end

"""
    Generates ω support of ρ for a given mesh type and model. 
"""
function generateω(mesh::LogMesh, model::PhysicalModel; offset=big(1e-15))
    # mesh parameters
    ω0, Nω, D = mesh.ω0, mesh.Nω, mesh.D
    Δ = model.Δ

    negfreqs = exp.(range(log(ω0), log(D - Δ), Int(Nω/2) - 1)) .+ Δ
    posfreqs = exp.(range(log(ω0),log(D - Δ), Int(Nω/2) - 1)) .+ Δ

    return vcat(-reverse(negfreqs), [-Δ - offset, Δ + offset], posfreqs)
end

########################
# DATA ANALYSIS
########################
"""
    Coefficients of the expansion of the D ∈ su(2) in the basis of Pauli matrices.
"""
function paulibasis(D::Matrix)
    σ1 = [0. 1.; 1. 0.]
    σ2 = [0. -im; im 0.]
    σ3 = [1. 0; 0. -1.]

    coeffs = zeros(ComplexF64, 4)
    coeffs[1] = tr(D)
    coeffs[2] = tr(σ1*D)
    coeffs[3] = tr(σ2*D)
    coeffs[4] = tr(σ3*D)

    return coeffs./2
end

"""
    Get the Pauli coefficients of a hybridization function for all ω.
"""
function getpaulicoeffs(hybri)
    Nω = size(hybri, 1)
    allcoeffs = zeros(4, Nω)
    for i in 1:Nω
        D = hybri[i, :, :]
        coeffs = paulibasis(D)
        allcoeffs[:, i] = real.(coeffs)
    end
    return real.(allcoeffs)
end

"""
    Get a tensor of values for a hybridization function ρ and an array of frequencies ωs. The shape of the tensor is (Nω, n, n), where Nω is the length of ωs and n is the number of bands. 
"""
function generateρs(ρ, ωs)
    ρs = [ρ(ω) for ω in ωs]
    n = size(first(ρs), 1)
    ρs = reduce(hcat, ρs)
    ρs = reshape(ρs, n, n, :)
    return permutedims(ρs, (3, 1, 2))
end

########################
# SAVING FILES
########################

"""
    Saves Wilson chains for each z in two different files: energies_i_Nz.txt and hoppings_i_Nz.txt. The default directory path is E_T_matrices/ (one can change this using the optional parameter foldername).
"""
function savechains(chains::Vector{WilsonChain};
                    foldername="E_T_matrices")
    mkpath(foldername)

    for chain in chains
        z, Nz = chain.z
        # on-site energies and hoppings
        E = chain.E
        T = chain.T
        # length of Wilson chain
        J = size(E, 1)

        suffix = "_$(z)_$(Nz).txt"

        # on-site energies
        energies_path = joinpath(foldername, "energies"*suffix)
        open(energies_path, "w") do f
            for i in 1:J
                e = round.(convert.(ComplexF64, E[i, :, :]), digits=13)
                writedlm(f, e)
                if i != J write(f, "\n") end
            end
        end

        # chain hoppings
        hoppings_path = joinpath(foldername, "hoppings"*suffix)
        open(hoppings_path, "w") do f
            for i in 1:J
                t = round.(convert.(ComplexF64, T[i, :, :]), digits=13)
                writedlm(f, t)
                if i != J write(f, "\n") end
            end
        end
    end
end

function loadchains()
    if isdir("E_T_matrices")
        files = readdir("E_T_matrices")
        Nz = Int(length(files)/2)
        chains = Vector{WilsonChain}(undef, Nz)
        for i in 1:Nz
            fn_e = "E_T_matrices/energies_$(i)_$(Nz).txt"
            fn_t = "E_T_matrices/hoppings_$(i)_$(Nz).txt"

            E = zeros(ComplexF64, 1, 1, 1)
            T = zeros(ComplexF64, 1, 1, 1)

            open(fn_e) do f
                lines = readlines(f)
                Nbands = length(split(lines[1], "\t"))
                J = Int((length(lines) + 1)/(Nbands + 1))
                E = zeros(ComplexF64, J, Nbands, Nbands)
                k = 1
                l = 1
                for (j, line) in enumerate(lines)
                    if line != ""
                        E[k, l, :] = parse.(ComplexF64, split(line, "\t"))
                        l += 1
                    else
                        k += 1
                        l = 1
                    end
                end
            end

            open(fn_t) do f
                lines = readlines(f)
                Nbands = length(split(lines[1], "\t"))
                J = Int((length(lines) + 1)/(Nbands + 1))
                T = zeros(ComplexF64, J, Nbands, Nbands)
                k = 1
                l = 1
                for (j, line) in enumerate(lines)
                    if line != ""
                        T[k, l, :] = parse.(ComplexF64, split(line, "\t"))
                        l += 1
                    else
                        k += 1
                        l = 1
                    end
                end
            end
            chains[i] = WilsonChain(E, T, (i, Nz))
        end
        return chains
    else
        error("Directory E_T_matrices does not exist in the working directory!")
    end

end