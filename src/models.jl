################################################
# PREDEFINED MODELS
################################################

abstract type PhysicalModel end

struct sWaveSC <: PhysicalModel
    # physical params
    # BCS gap
    Δ
    # skewness of the density of states
    κ
    # hybridization strength
    Γ
end

sWaveSC() = sWaveSC(0.1, 0., 0.5)
sWaveSC(Δ) = sWaveSC(Δ, 0., 0.5)

struct Flat{T<:Real} <: PhysicalModel
    Δ::T
    Γ::T
end

Flat() = Flat(0., 0.5/pi)

# Hybridization function of the double SC lead AIM 
# in the scalar representation
# Taken from Zalom's paper arXiv:2307.07479 Eq. (19)
# Compared to sWaveSc it has an additional parameter ϕ
# which is the difference of SC phases in both leads
struct DoubleLeadSC <: PhysicalModel
    # BCS gap
    Δ::Float64
    # hybridization strength
    Γ::Float64
    # phase difference
    ϕ::Float64
end

DoubleLeadSC() = DoubleLeadSC(0.01, 0.5, 0.)
DoubleLeadSC(Δ, ϕ) = DoubleLeadSC(Δ, 0.5, ϕ)

################################################
# HYBRIDIZATIONS FOR PREDEFINED MODELS
################################################

function hybridization(ω::Float64, model::sWaveSC; minvalue=0.)
    Δ, κ, Γ = model.Δ, model.κ, model.Γ
    if abs(ω) < Δ
        return [minvalue minvalue; minvalue minvalue]
    else
        ξ = sqrt(ω^2 - Δ^2)
        τ₀ = [1. 0.; 0. 1.]
        τ₁ = [0. 1.; 1. 0.]
        τ₃ = [1. 0.; 0. -1.]
        return (Γ/pi) * (1/ξ) * (abs(ω)*τ₀ + κ*ω*ξ*τ₃ - Δ*sign(ω)*τ₁)
    end
end

function hybridization(ω::BigFloat, model::sWaveSC; minvalue=big"0.")
    Δ, κ, Γ = model.Δ, model.κ, model.Γ
    if abs(ω) < Δ
        return [minvalue minvalue; minvalue minvalue]
    else
        ξ = sqrt(ω^2 - Δ^2)
        τ₀ = [big"1." big"0."; big"0." big"1."]
        τ₁ = [big"0." big"1."; big"1." big"0."]
        τ₃ = [big"1." big"0."; big"0." big"-1."]
        return (Γ/pi) * (1/ξ) * (abs(ω)*τ₀ + κ*ω*ξ*τ₃ - Δ*sign(ω)*τ₁)
    end
end

struct Bethe
    t
end

function hybridization(ω, model::Bethe; minvalue=big"0.")
    t = model.t
    if abs(ω) < 2t
        return (1/(2t^2))*sqrt(4t^2 - ω^2)
    else
        return minvalue
    end
end