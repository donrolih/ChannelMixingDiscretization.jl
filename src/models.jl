########################
# PREDEFINED MODELS
########################

struct sWaveSC{T<:Real}
    # physical params
    # BCS gap
    Δ::T
    # impurity level
    εf::T
    # skewness of the density of states
    κ::T
    # Hubbard repulsion
    U::T
    # hybridization strength
    Γ::T
end

sWaveSC() = sWaveSC(0.1, 0., 0., 0., 0.5/pi)

struct Flat{T<:Real}
    Δ::T
    Γ::T
end

Flat() = Flat(0., 0.5/pi)