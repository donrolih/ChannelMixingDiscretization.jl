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
using QuadGK
using DataInterpolations

# precision when using BigFloats
setprecision(256)

# from meshes.jl
export LogMesh, LinMesh, DiscretizationParams
# from models.jl
export sWaveSC, Flat, DoubleLeadSC, hybri
# from chain.jl
export WilsonChain, maptochains
# from star.jl
export StarHamiltonian, discmodel
# from utils.jl
export twistingparameters, paulibasis, generateω, getpaulicoeffs, generateρs, savechains, loadchains, discretize
# from reconstruct.jl
export reconstructhybri, spectralfunction, broadenaverage

include("meshes.jl")
include("models.jl")
include("star.jl")
include("tridiagonalization.jl")
include("chain.jl")
include("utils.jl")
include("reconstruct.jl")
include("integrators.jl")

end