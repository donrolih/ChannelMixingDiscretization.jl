function discretegrid(x₀, Λ, N, z)
    # Initialize an array of zeros with length N
    xs = zeros(typeof(x₀), N)
    
    # Fill the array with values based on the geometric progression
    for i in 1:N
        x = i + z
        if x ≤ 2
            xs[i] = x₀
        else
            xs[i] = x₀ * Λ^(2 - x)
        end
    end
    
    # Return the resulting array
    return xs
end

function getweights(f::DataInterpolations.LinearInterpolation, ωs)
    N = length(ωs)
    T = typeof(first(ωs))
    posweights = zeros(T, N - 1)
    negweights = zeros(T, N - 1)
    for n in 1:(N-1)
        # we are iterating from the right side of interval: x[i+1] < x[i]
        a = ωs[n+1]
        b = ωs[n]
        posweights[n] = DataInterpolations.integral(f, a, b)
        negweights[n] = DataInterpolations.integral(f, -b, -a)
    end
    return negweights, posweights
end

function getenergies(ωs)
    N = length(ωs)
    E = zeros(typeof(first(ωs)), N-1)
    for n in 1:(N-1)
        E[n] = (ωs[n] + ωs[n+1]) / 2
    end
    return E
end

function bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N, z)
    Δ_interp = DataInterpolations.LinearInterpolation(Δ, input_ωs; extrapolate=true)
    Δoff_interp = DataInterpolations.LinearInterpolation(Δoff, input_ωs; extrapolate=true)
    
    D = input_ωs[end]
    ωs = discretegrid(D, Λ, N, z)

    w_neg, w_pos = getweights(Δ_interp, ωs)
    wbar_neg, wbar_pos = getweights(Δoff_interp, ωs)

    E = getenergies(ωs)

    δ_pos = (wbar_pos ./ w_pos) .* E
    δ_neg = -(wbar_neg ./ w_neg) .* E

    ξ_pos = sqrt.(1 .- (wbar_pos .^ 2) ./ (w_pos .^ 2)) .* E
    ξ_neg = -sqrt.(1 .- (wbar_neg .^ 2) ./ (w_neg .^ 2)) .* E

    γ_pos = sqrt.(w_pos)
    γ_neg = sqrt.(w_neg)

    return δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg
end


function bauer_discretization_chain(input_ωs, Δ, Δoff, Λ, N, Nmax, z)
    δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg = bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N, z)

    Nstar = length(δ_pos)
    @assert Nmax < N "Nmax must be smaller than N"

    # we start by defining the localised state's coupling to the impurity
    # called "theta" in NRG Ljubljana 
    β_min = sqrt(sum(γ_pos .^ 2 + γ_neg .^ 2)) # Eq. (29)
    T = typeof(β_min)
    @show T

    # diagonal hoppings
    βs = zeros(T, Nmax)

    # on-site energies
    εs = zeros(T, Nmax)
    # SC pairing
    Δs = zeros(T, Nmax)

    U = Dict{Int, Array{T, 2}}()
    U[1] = zeros(T, Nmax, Nstar)
    U[-1] = zeros(T,Nmax, Nstar)
    
    V = Dict{Int, Array{T, 2}}()
    V[1] = zeros(T, Nmax, Nstar)
    V[-1] = zeros(T, Nmax, Nstar)
    
    # Initialize the values of U and V (V is zero)
    U[1][1, :] = γ_pos ./ β_min
    U[-1][1, :] = γ_neg ./ β_min

    # @show U[1][1, :]
    # @show U[-1][1, :]

    # Eq. (33)
    ε_pos = ξ_pos .* (U[1][1, :] .^ 2 .- V[1][1, :] .^ 2) + 2δ_pos .* U[1][1, :] .* V[1][1, :]
    ε_neg = ξ_neg .* (U[-1][1, :] .^ 2 .- V[-1][1, :] .^ 2) + 2δ_neg .* U[-1][1, :] .* V[-1][1, :]
    εs[1] = sum(ε_pos) + sum(ε_neg)

    # Eq. (34)
    Δ_pos = δ_pos .* (U[1][1, :] .^ 2 .- V[1][1, :] .^ 2) - 2ξ_pos .* U[1][1, :] .* V[1][1, :]
    Δ_neg = δ_neg .* (U[-1][1, :] .^ 2 .- V[-1][1, :] .^ 2) - 2ξ_neg .* U[-1][1, :] .* V[-1][1, :]
    Δs[1] = sum(Δ_pos) + sum(Δ_neg)

    # Eq. (35): for n = 1 (n=0 in Bauer notation) there is no additional β^2 term on the RHS
    β_pos = (ξ_pos .^ 2 .+ δ_pos .^ 2) .* (U[1][1, :] .^ 2 .+ V[1][1, :] .^ 2)
    β_neg = (ξ_neg .^ 2 .+ δ_neg .^ 2) .* (U[-1][1, :] .^ 2 .+ V[-1][1, :] .^ 2)
    βs[1] = sum(β_pos) + sum(β_neg) - εs[1]^2 - Δs[1]^2
    @show convert.(Float64, βs[1])
    βs[1] = sqrt(βs[1])


    for n in 2:Nmax
        @info "Current index: " n
        if n==2
            U[1][n, :] = (ξ_pos .- εs[n-1]) .* U[1][n-1, :] + (δ_pos .+ Δs[n-1]) .* V[1][n-1, :]
            U[1][n, :] =  U[1][n, :] ./ βs[n-1]
            U[-1][n, :] = (ξ_neg .- εs[n-1]) .* U[-1][n-1, :] + (δ_neg .+ Δs[n-1]) .* V[-1][n-1, :]
            U[-1][n, :] = U[-1][n, :] ./ βs[n-1]

            V[1][n, :] = (δ_pos .- Δs[n-1]) .* U[1][n-1, :] - (ξ_pos .+ εs[n-1]) .* V[1][n-1, :]
            V[1][n, :] = V[1][n, :] ./ βs[n-1]
            V[-1][n, :] = (δ_neg .- Δs[n-1]) .* U[-1][n-1, :] - (ξ_neg .+ εs[n-1]) .* V[-1][n-1, :]
            V[-1][n, :] =  V[-1][n, :] ./ βs[n-1]
        else
            U[1][n, :] = (ξ_pos .- εs[n-1]) .* U[1][n-1, :] + (δ_pos .+ Δs[n-1]) .* V[1][n-1, :] - βs[n-2] .* U[1][n-2, :]
            U[1][n, :] =  U[1][n, :] ./ βs[n-1]
            U[-1][n, :] = (ξ_neg .- εs[n-1]) .* U[-1][n-1, :] + (δ_neg .+ Δs[n-1]) .* V[-1][n-1, :] - βs[n-2] .* U[-1][n-2, :]
            U[-1][n, :] = U[-1][n, :] ./ βs[n-1]

            V[1][n, :] = (δ_pos .- Δs[n-1]) .* U[1][n-1, :] - (ξ_pos .+ εs[n-1]) .* V[1][n-1, :] - βs[n-2] .* V[1][n-2, :]
            V[1][n, :] = V[1][n, :] ./ βs[n-1]
            V[-1][n, :] = (δ_neg .- Δs[n-1]) .* U[-1][n-1, :] - (ξ_neg .+ εs[n-1]) .* V[-1][n-1, :] - βs[n-2] .* V[-1][n-2, :]
            V[-1][n, :] =  V[-1][n, :] ./ βs[n-1]
        end

        ε_pos = ξ_pos .* (U[1][n, :] .^ 2 .- V[1][n, :] .^ 2) + 2δ_pos .* U[1][n, :] .* V[1][n, :]
        ε_neg = ξ_neg .* (U[-1][n, :] .^ 2 .- V[-1][n, :] .^ 2) + 2δ_neg .* U[-1][n, :] .* V[-1][n, :]
        εs[n] = sum(ε_pos) + sum(ε_neg)

        Δ_pos = δ_pos .* (U[1][n, :] .^ 2 .- V[1][n, :] .^ 2) - 2ξ_pos .* U[1][n, :] .* V[1][n, :]
        Δ_neg = δ_neg .* (U[-1][n, :] .^ 2 .- V[-1][n, :] .^ 2) - 2ξ_neg .* U[-1][n, :] .* V[-1][n, :]
        Δs[n] = sum(Δ_pos) + sum(Δ_neg)

        β_pos = (ξ_pos .^ 2 .+ δ_pos .^ 2) .* (U[1][n, :] .^ 2 .+ V[1][n, :] .^ 2)
        β_neg = (ξ_neg .^ 2 .+ δ_neg .^ 2) .* (U[-1][n, :] .^ 2 .+ V[-1][n, :] .^ 2)

        βs[n] = sum(β_pos) + sum(β_neg) - εs[n]^2 - Δs[n]^2 - βs[n-1]^2
        βs[n] = sqrt(βs[n])
    end

    return U, V, εs, Δs, βs, β_min
end

function bauer_discretization_chain_matrix(input_ωs, Δ, Δoff, Λ, N, Nmax, Nz)
    chains = Vector{WilsonChain}(undef, Nz)
    for i in 1:Nz
        z = i / Nz
        @info "Mapping to a Wilson chain for twisting parameter z = $(z) ..."
        δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg = bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N, z)
        Nstar = length(δ_pos)
        @assert Nmax < N "Nmax must be smaller than N"
        T = typeof(first(δ_pos))
        X_plus = zeros(T, 2, 2)
        X_minus = zeros(T, 2, 2)

        X_plus[1, 1] = ξ_pos[1]
        X_plus[1, 2] = -δ_pos[1]
        X_plus[2, 1] = -δ_pos[1]
        X_plus[2, 2] = -ξ_pos[1]

        X_minus[1, 1] = ξ_neg[1]
        X_minus[1, 2] = -δ_neg[1]
        X_minus[2, 1] = -δ_neg[1]
        X_minus[2, 2] = -ξ_neg[1]
        # X_minus = -X_minus

        Γ_plus = zeros(T, 2, 2)
        Γ_minus = zeros(T, 2, 2)

        Γ_plus[1, 1] = γ_pos[1]
        Γ_plus[2, 2] = γ_pos[1]

        Γ_minus[1, 1] = γ_neg[1]
        Γ_minus[2, 2] = γ_neg[1]

        ζ = zeros(T, 2, 2)
        ζ += sum((γ_pos .^ 2 + γ_neg .^ 2)) * LinearAlgebra.I
        @info "ζ = $ζ"
        U = inv(sqrt(ζ)) * Γ_plus
        V = inv(sqrt(ζ)) * Γ_minus

        E0 = U * X_plus * U' + V * X_minus * V'


        T0_sqr = (U * X_plus' - E0' * U) * (X_plus * U' - U' * E0) + (V * X_minus' - E0' * V) * (X_minus * V' - V' * E0)

        T0 = sqrt(T0_sqr)

        Ts = zeros(T, (Nmax + 1, 2, 2))
        Es = zeros(T, (Nmax, 2, 2))

        Ts[1, :, :] = ζ
        Ts[2, :, :] = T0
        
        Es[1, :, :] = E0

        for i=2:Nmax
            U = (X_plus * U' - U' * E0) * inv(T0)
            V = (X_minus * V' - V' * E0) * inv(T0)

            X_plus[1, 1] = ξ_pos[i]
            X_plus[1, 2] = -δ_pos[i]
            X_plus[2, 1] = -δ_pos[i]
            X_plus[2, 2] = -ξ_pos[i]

            X_minus[1, 1] = ξ_neg[i]
            X_minus[1, 2] = -δ_neg[i]
            X_minus[2, 1] = -δ_neg[i]
            X_minus[2, 2] = -ξ_neg[i]

            E1 = U * X_plus * U' + V * X_minus * V'
            T1_sqr = (U * X_plus' - E1' * U) * (X_plus * U' - U' * E1) + (V * X_minus' - E1' * V) * (X_minus * V' - V' * E1)
            T1 = sqrt(T1_sqr)

            Ts[i + 1, :, :] = T1
            Es[i, :, :] = E1

            E0 = E1
            T0 = T1
        end

        chains[i] = WilsonChain(abs.(convert.(ComplexF64, Es)), convert.(ComplexF64, Ts), (1, 1))
    end
    return chains
end

function bauer_discretization_chain_matrix_second(input_ωs, Δ, Δoff, Λ, N, Nmax, Nz)
    chains = Vector{WilsonChain}(undef, Nz)
    for i in 1:Nz
        z = i / Nz
        @info "Mapping to a Wilson chain for twisting parameter z = $(z) ..."
        δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg = bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N, z)
        Nstar = length(δ_pos)
        @assert Nmax < N "Nmax must be smaller than N"
        T = typeof(first(γ_pos))

        X_plus = zeros(T, Nstar, 2, 2)
        X_minus = zeros(T, Nstar, 2, 2)

        X_plus[:, 1, 1] = ξ_pos
        X_plus[:, 1, 2] = -δ_pos
        X_plus[:, 2, 1] = -δ_pos
        X_plus[:, 2, 2] = -ξ_pos

        X_minus[:, 1, 1] = ξ_neg
        X_minus[:, 1, 2] = -δ_neg
        X_minus[:, 2, 1] = -δ_neg
        X_minus[:, 2, 2] = -ξ_neg

        U = zeros(T, Nmax, Nstar, 2, 2)
        V = zeros(T, Nmax, Nstar, 2, 2)

        Γ_plus = zeros(T, Nstar, 2, 2)
        Γ_minus = zeros(T, Nstar, 2, 2)

        Γ_plus[:, 1, 1] = γ_pos
        Γ_plus[:, 2, 2] = γ_pos

        Γ_minus[:, 1, 1] = γ_neg
        Γ_minus[:, 2, 2] = γ_neg

        Es = zeros(T, Nmax, 2, 2)
        Ts = zeros(T, Nmax, 2, 2)

        ζ = zeros(2, 2)
        for j in 1:Nstar
            ζ += Γ_plus[j, :, :] * Γ_plus[j, :, :] + Γ_minus[j, :, :] * Γ_minus[j, :, :]
        end
        @show convert.(Float64, sqrt(ζ))
        @show convert.(Float64, inv(sqrt(ζ)))

        T0_sqr = zeros(T, 2, 2)
        for j in 1:Nstar
            U[1, j, :, :] = inv(sqrt(ζ)) * Γ_plus[j, :, :]
            V[1, j, :, :] = inv(sqrt(ζ)) * Γ_minus[j, :, :]
            Es[1, :, :] += U[1, j, :, :] * X_plus[j, :, :] * U[1, j, :, :]' + V[1, j, :, :] * X_minus[j, :, :] * V[1, j, :, :]'
            M0 = U[1, j, :, :] * X_plus[j, :, :]' - Es[1, :, :]' * U[1, j, :, :]
            N0 = V[1, j, :, :] * X_minus[j, :, :]' - Es[1, :, :]' * V[1, j, :, :]
            T0_sqr += M0 * M0' + N0 * N0'
        end
        @show convert.(Float64, U)
        @show convert.(Float64, Es[1, :, :])
        @show convert.(Float64, T0_sqr)
        @show convert.(Float64, sqrt(T0_sqr))
        Ts[1, :, :] = sqrt(T0_sqr)

        for j in 1:Nstar
            U[2, j, :, :] = ((X_plus[j, :, :] * U[1, j, :, :]' - U[1, j, :, :]' * Es[1, :, :]) * inv(Ts[1, :, :]'))'
            V[2, j, :, :] = ((X_minus[j, :, :] * V[1, j, :, :]' - V[1, j, :, :]' * Es[1, :, :]) * inv(Ts[1, :, :]'))'
        end

        for i=2:Nmax
            @info "Current index: " i
            for j in 1:Nstar
                Es[i, :, :] += U[i, j, :, :] * X_plus[j, :, :] * U[i, j, :, :]' + V[i, j, :, :] * X_minus[j, :, :] * V[i, j, :, :]'
            end

            T_sqr = zeros(T, 2, 2)
            for j in 1:Nstar
                Mj = U[i, j, :, :] * X_plus[j, :, :]' - Es[i, :, :]' * U[i, j, :, :] - Ts[i-1, :, :]' * U[i-1, j, :, :]
                Nj = V[i, j, :, :] * X_minus[j, :, :]' - Es[i, :, :]' * V[i, j, :, :] - Ts[i-1, :, :]' * V[i-1, j, :, :]
                T_sqr += Mj * Mj' + Nj * Nj'
            end
            Ts[i, :, :] = sqrt(T_sqr)
            
            if i == Nmax
                break
            else
                for j in 1:Nstar
                    U[i+1, j, :, :] = ((X_plus[j, :, :] * U[i, j, :, :]' - U[i, j, :, :]' * Es[i, :, :] - U[i-1, j, :, :]' * Ts[i - 1, :, :]) * inv(Ts[i, :, :]'))'
                    V[i+1, j, :, :] = ((X_minus[j, :, :] * V[i, j, :, :]' - V[i, j, :, :]' * Es[i, :, :] - V[i-1, j, :, :]' * Ts[i - 1, :, :]) * inv(Ts[i, :, :]'))'
                end
            end
        end

        chains[i] = WilsonChain(abs.(convert.(ComplexF64, Es)), convert.(ComplexF64, Ts), (1, 1))
    end
    return chains
end