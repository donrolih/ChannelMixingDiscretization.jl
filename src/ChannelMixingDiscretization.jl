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
export logmesh, LinMesh, DiscretizationParams
# from models.jl
export sWaveSC, Flat, DoubleLeadSC, Bethe, hybridization
# from chain.jl
export WilsonChain, maptochains
# from starrepresentation.jl
export StarHamiltonian
# from wrapper.jl
export discretize
# from utils.jl
export twistingparameters, paulibasis, generateω, getpaulicoeffs, generateρs, savechains, loadchains
# from reconstruct.jl
export reconstructhybri, spectralfunction, broadenaverage

include("meshes.jl")
include("models.jl")
include("starrepresentation.jl")
include("tridiagonalization.jl")
include("chain.jl")
include("wrapper.jl")
include("utils.jl")
include("reconstruct.jl")
include("integrators.jl")

end