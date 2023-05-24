module ChannelMixingDiscretization

# Write your package code here.
# Dependencies
using LinearAlgebra
using GenericLinearAlgebra
using BlockDiagonals
using BlockBandedMatrices
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
export twistingparameters, paulibasis, generateω, getpaulicoeffs, generateρs, savechains, loadchains
# from reconstruct.jl
export reconstructhybri, spectralfunction

include("meshes.jl")
include("models.jl")
include("star.jl")
include("tridiagonalization.jl")
include("chain.jl")
include("utils.jl")
include("reconstruct.jl")

end