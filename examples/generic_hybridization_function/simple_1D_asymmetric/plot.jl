##
using ChannelMixingDiscretization
using PyPlot, LaTeXStrings
plt.style.use(["don-custom"])

# function to discretize
ρ(ω; κ=1.) = (1/2)*(1 + κ*ω)

# mesh parameters

ω0 = 1e-4
Nω = 5000
D = 1.

mesh = LogMesh(ω0, Nω, D)
model = Flat()

# discretization parameters

Nz = 32
J = 60
zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))

params = DiscretizationParams(2., zs, J, 500000, 1e-5);

# mapping to star model
starH = discmodel(ρ, mesh, model, params);

# mapping to Wilson chain
chains = maptochains(starH; m=J);

ωs = generateω(mesh, model; offset=1e-4)
reconstructed_chain = reshape(reconstructhybri(chains, ωs; smear=0.7), Nω);

##
fig, ax = plt.subplots(nrows=2, figsize=(12, 5), sharex=true)

ax[1].plot(ωs, reconstructed_chain, ".")
ax[1].plot(ωs, ρ.(ωs))

dif = abs.(reconstructed_chain .- ρ.(ωs))
ax[2].plot(ωs, dif)

ax[2].set_yscale("log")

fig.tight_layout()

fig.savefig("wilson_chain_mapping.pdf")
