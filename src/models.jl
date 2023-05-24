################################################
# PREDEFINED MODELS
################################################

abstract type PhysicalModel end

struct sWaveSC <: PhysicalModel
    # physical params
    # BCS gap
    Δ::Float64
    # impurity level
    εf::Float64
    # skewness of the density of states
    κ::Float64
    # Hubbard repulsion
    U::Float64
    # hybridization strength
    Γ::Float64
end

sWaveSC() = sWaveSC(0.1, 0., 0., 0., 0.5/pi)
sWaveSC(Δ) = sWaveSC(Δ, 0., 0., 0., 0.5/pi)

struct Flat{T<:Real} <: PhysicalModel
    Δ::T
    Γ::T
end

Flat() = Flat(0., 0.5/pi)

################################################
# PREDEFINED HYBRIDIZATIONS
################################################

function hybri(model::sWaveSC, mesh::Mesh; η=1e-10)
    function hybri(ω)
        σ0 = [1. 0.; 0. 1.]
        σ1 = [0. 1.; 1. 0.]
        σ3 = [1. 0.; 0. -1.]

        Δ = model.Δ
        Γ = model.Γ
        κ = model.κ
        
        D = mesh.D
        z = ω + η*im
        ξ = sqrt(Δ^2 - z^2)

        I0 = -2atan(D/ξ)/ξ
        I2 = -2D + 2ξ*atan(D/ξ)
        Σ = (Γ/pi)*(z*I0*σ0 - I0*Δ*σ1 + κ*I2*σ3)
        return (im/(2pi))*(Σ - Σ')
    end
end

function hybri(model::Flat, mesh::Mesh)
    hybri(ω) = model.Γ
end