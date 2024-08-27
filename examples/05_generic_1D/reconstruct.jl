using ChannelMixingDiscretization
using DelimitedFiles
using DataInterpolations
using Integrals
using Plots

# load the Wilson chains from E_T_matrices folder
chains = loadchains()

# reconstruct the original hybridization function
lb = -5
ub = 5
Nω = 5000

# broadening parameter
η = 0.15

ωs, hyb = broaden1D(chains, lb, ub, Nω; η=η)

data = readdlm("DOS.dat")


problem = SampledIntegralProblem(-data[:, 2], data[:, 1])
method = TrapezoidalRule()
sol = solve(problem, method)

norm_dos = sol.u
hybridization = LinearInterpolation(-data[:, 2] ./ norm_dos, data[:, 1])

p = plot()
plot!(p, size=(1000, 750))
plot!(p, ωs, real(hyb), label="reconstruted")
plot!(p, ωs, hybridization(ωs), label="input")
savefig(p, "reconstructed.pdf")