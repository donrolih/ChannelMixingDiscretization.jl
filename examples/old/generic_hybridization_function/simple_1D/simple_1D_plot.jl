using ChannelMixingDiscretization
##
# offset is important (ω0!!! if it is too small, you get divergences at zero!)
mesh = LogMesh(1e-4, 5000, 1.)
model = Flat()

Nz = 32
J = 60
zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
params = DiscretizationParams(1.8, zs, J, 500000, 1e-5);

# hybridization function to map
ρ = hybri(model, mesh);
# mapping to star model
starH = discmodel(ρ, mesh, model, params);
##
ωs = generateω(mesh, model; offset=1e-4)
ρs = generateρs(ρ, ωs)
reconstructed_star = reconstructhybri(starH, ωs; smear=0.8);
##
using PyPlot, LaTeXStrings
plt.style.use("don-custom")
##
fig, ax = plt.subplots(nrows=2, sharex=true)
ρs = reshape(ρs, 5000)
ax[1].plot(ωs, ρs)
ax[1].plot(ωs, reconstructed_star, ".")

dif = abs.(ρs .- reconstructed_star)
ax[2].plot(ωs, dif)
ax[2].set_yscale("log")

ax[1].set_ylabel(L"$\mathcal{D}(\omega)$")
ax[2].set_xlabel(L"$\omega$")
ax[2].set_ylabel("error")

plt.tight_layout()

Nz = length(params.z)
Λ = params.Λ
title = L"Mapping to a star Hamiltonian: $N_z = %$(Nz)$, $\Lambda = %$(Λ)$"
fig.suptitle(y=1.01, title)
display(fig)
##
# mapping to Wilson chains
chains = maptochains(starH; m=J);
reconstructed_chain = reshape(reconstructhybri(chains, ωs; smear=0.7), 5000);
##
fig, ax = plt.subplots(nrows=2, sharex=true)
ρs = reshape(ρs, 5000)
ax[1].plot(ωs, ρs)
ax[1].plot(ωs, reconstructed_chain, ".")

dif = abs.(ρs .- reconstructed_chain)
ax[2].plot(ωs, dif)
ax[2].set_yscale("log")

ax[1].set_ylabel(L"$\mathcal{D}(\omega)$")
ax[2].set_xlabel(L"$\omega$")
ax[2].set_ylabel("error")

plt.tight_layout()

Nz = length(params.z)
Λ = params.Λ
title = L"Mapping to a Wilson chain: $N_z = %$(Nz)$, $\Lambda = %$(Λ)$"
fig.suptitle(y=1.01, title)
display(fig)
# fig.savefig("wilson_chain_mapping.pdf")
##
# spectral function with diagonalisation
data = spectralfunction(chains)
energies = data[:, 1]
spectral = data[:, 2];

##
fig, ax = plt.subplots()

ax.plot(energies, spectral, ".")

display(fig)