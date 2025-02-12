function discretegrid(x₀, Λ, N)
    # Initialize an array of zeros with length N
    xs = zeros(N)
    
    # Fill the array with values based on the geometric progression
    for i in 1:N
        xs[i] = x₀ * Λ^(-i+1)
    end
    
    # Return the resulting array
    return xs
end

function getweights(f::DataInterpolations.LinearInterpolation, ωs)
    N = length(ωs)
    posweights = zeros(N - 1)
    negweights = zeros(N - 1)
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
    E = zeros(N-1)
    for n in 1:(N-1)
        E[n] = (ωs[n] + ωs[n+1]) / 2
    end
    return E
end

function bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N)
    Δ_interp = DataInterpolations.LinearInterpolation(Δ, input_ωs; extrapolate=true)
    Δoff_interp = DataInterpolations.LinearInterpolation(Δoff, input_ωs; extrapolate=true)
    
    D = input_ωs[end]
    ωs = discretegrid(D, Λ, N)

    w_neg, w_pos = getweights(Δ_interp, ωs)
    wbar_neg, wbar_pos = getweights(Δoff_interp, ωs)

    E = getenergies(ωs)

    δ_pos = (wbar_pos ./ w_pos) .* E
    δ_neg = -(wbar_neg ./ w_neg) .* E

    x_pos = (wbar_pos .^ 2) ./ (w_pos .^ 2)
    x_neg = (wbar_neg .^ 2) ./ (w_neg .^ 2)
    
    @info x_pos
    @info x_neg

    ξ_pos = sqrt.(1 .- (wbar_pos .^ 2) ./ (w_pos .^ 2)) .* E
    ξ_neg = -sqrt.(1 .- (wbar_neg .^ 2) ./ (w_neg .^ 2)) .* E

    γ_pos = sqrt.(w_pos)
    γ_neg = sqrt.(w_neg)

    return δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg
end


function bauer_discretization_chain(input_ωs, Δ, Δoff, Λ, N, Nmax)
    δ_pos, δ_neg, ξ_pos, ξ_neg, γ_pos, γ_neg = bauer_discretization_star(input_ωs, Δ, Δoff, Λ, N)

    N = length(δ_pos)
    @assert Nmax < N "Nmax must be smaller than N"

    # we start by defining the localised state's coupling to the impurity
    # called "theta" in NRG Ljubljana 
    β_min = sqrt(sum(γ_pos .^ 2 + γ_neg .^ 2)) # Eq. (29)

    # diagonal hoppings
    βs = zeros(Nmax)

    # on-site energies
    εs = zeros(Nmax)
    # SC pairing
    Δs = zeros(Nmax)

    U = Dict{Int, Array{T, 2}}()
    V = Dict{Int, Array{T, 2}}()

    # Eq. (33)
    ε_pos = ξ_pos .^ 2 .* (U[1][1, :] .^ 2 .- V[1][1, :] .^ 2) + 2δ_neg .* U[-1][1, :] .* V[-1][1, :]
    ε_neg = ξ_neg .^ 2 .* (U[-1][1, :] .^ 2 .- V[-1][1, :] .^ 2) + 2δ_neg .* U[-1][1, :] .* V[-1][1, :]
    εs[1] = sum(ε_pos) + sum(ε_neg)

    # Eq. (34)
    Δ_neg = δ_neg .* (U[-1][1, :] .^ 2 .- V[-1][1, :] .^ 2) + 2ξ_neg .* U[-1][1, :] .* V[-1][1, :]
    Δ_pos = δ_pos .* (U[1][1, :] .^ 2 .- V[1][1, :] .^ 2) + 2ξ_pos .* U[1][1, :] .* V[1][1, :]
    Δs[1] = sum(Δ_pos) + sum(Δ_neg)

    # Eq. (35)
    
end