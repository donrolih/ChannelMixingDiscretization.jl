module ChannelMixingDiscretization

# Write your package code here.
# Dependencies
using LinearAlgebra
using BlockDiagonals
using NumericalIntegration
using Interpolations
using DelimitedFiles, JLD2

# from meshes.jl
export LogMesh, DiscretizationParams
# from models.jl
export sWaveSC, Flat, hybri
# from chain.jl
export WilsonChain, maptochains
# from star.jl
export StarHamiltonian, discmodel
# from utils.jl
export twistingparameters, paulibasis, generateω, getpaulicoeffs, generateρs
# from reconstruct.jl
export reconstructhybri
export interpolate, int

function interpolate()
    R = 1:40
    ωs = 1:40
    return linear_interpolation(R, ωs, extrapolation_bc=Line())
end

function int()
    ωs = 1:40
    weight = 1:40
    return cumul_integrate(ωs, weight)
end

include("meshes.jl")
include("models.jl")
include("star.jl")
include("tridiagonalization.jl")
include("chain.jl")
include("utils.jl")
include("reconstruct.jl")

end