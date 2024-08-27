using ChannelMixingDiscretization
using Plots

# load the Wilson chains from E_T_matrices folder
chains = loadchains()

# reconstruct the original hybridization function
lb = -1.5
ub = 1.5
Nω = 5000

# broadening parameter
η = 0.15

ωs, hyb_diag, hyb_offdiag = broadenSC(chains, lb, ub, Nω; η=η)

p = plot()
plot!(p, size=(1000, 750))
plot!(p, ωs, real(hyb_diag), label="diag")
plot!(p, ωs, real(hyb_offdiag), label="diag")
savefig(p, "reconstructed.pdf")