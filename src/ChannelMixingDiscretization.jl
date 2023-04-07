module ChannelMixingDiscretization

# Write your package code here.
# Dependencies
using LinearAlgebra
using BlockDiagonals
using NumericalIntegration
using Interpolations

# from meshes.jl
export LogMesh, DiscretizationParams
# from models.jl
export sWaveSC, Flat
# from chain.jl
export WilsonChain, maptochains
# from star.jl
export StarModel, maptostar
# from utils.jl
export twistingparameters, paulibasis, generateω
# from reconstruct.jl
export reconstructhybri

include("meshes.jl")
include("models.jl")
include("star.jl")
include("chain.jl")
include("discretization.jl")
include("utils.jl")
include("reconstruct.jl")

end