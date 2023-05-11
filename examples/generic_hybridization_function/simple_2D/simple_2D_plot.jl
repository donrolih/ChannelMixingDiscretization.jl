using ChannelMixingDiscretization
##
# offset is important (ω0!!! if it is too small, tou get divergences at zero!)
mesh = LogMesh(1e-4, 1000, 1.)
model = sWaveSC(0., 0., 0., 0., 0.5/pi)
Nz = 64
zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
J = 35
params = DiscretizationParams(2., zs, J, 500000, 1e-5)
##
function ρ(ω::Real)
    σ0 = [1. 0; 0. 1.]
    σ1 = [0. 1.; 1. 0.]
    σ2 = [0. -im ; im 0.]
    return σ0 + 0.4*σ1 + 0.5*(0.1 + ω^2)*σ2
end

starH = discmodel(ρ, mesh, model, params);
##
ωs = generateω(mesh, model; offset=1e-4)
reconstructed_star = reconstructhybri(starH, ωs; smear=0.7);
##
ρs = generateρs(ρ, ωs)
pauli_coeffs = getpaulicoeffs(ρs)
pauli_coeffs_re = getpaulicoeffs(reconstructed_star);
##
using PyPlot, LaTeXStrings
plt.style.use("don-custom")
##
fig, ax = plt.subplots(nrows=2, sharex=true)

labels = [L"$d_0$", L"$d_x$", L"$d_y$", L"$d_z$"]
for i in 1:4
    ax[1].plot(ωs, pauli_coeffs[i, :], label=labels[i], color="C$(i-1)")
    ax[1].plot(ωs, pauli_coeffs_re[i, :], ".", label=labels[i], color="C$(i-1)")
    dif = abs.(pauli_coeffs_re[i, :] .- pauli_coeffs[i, :])
    ax[2].plot(ωs, dif, label=labels[i], color="C$(i-1)")
end
ax[2].axhline(y=1e-2, ls="--", c="gray")

ax[2].set_yscale("log")

# labels
ax[1].set_ylabel(L"$\mathcal{D}(\omega)$")
ax[2].set_xlabel(L"$\omega$")
ax[2].set_ylabel("error")

ax[2].legend(ncol=2)
plt.tight_layout()
Nz = length(params.z)
Λ = params.Λ
title = L"Mapping to a star Hamiltonian: $N_z = %$(Nz)$, $\Lambda = %$(Λ)$"
fig.suptitle(y=1.01, title)
display(fig)

##
# from here on: mapping to the Wilson chain
chains = maptochains(starH; m=J);
##
reconstructed_chain = reconstructhybri(chains, ωs; smear=0.7)
pauli_coeffs_chain = getpaulicoeffs(reconstructed_chain);
##
fig, ax = plt.subplots(nrows=2, sharex=true)

labels = [L"$d_0$", L"$d_x$", L"$d_y$", L"$d_z$"]
for i in 1:4
    ax[1].plot(ωs, pauli_coeffs[i, :], label=labels[i], color="C$(i-1)")
    ax[1].plot(ωs, pauli_coeffs_chain[i, :], ".", label=labels[i], color="C$(i-1)")
    dif = abs.(pauli_coeffs_chain[i, :] .- pauli_coeffs[i, :])
    ax[2].plot(ωs, dif, label=labels[i], color="C$(i-1)")
end

ax[1].set_ylim(-0.1, 1.2)
# ax[1].set_xlim(-1e-4, 1e-4)

ax[2].axhline(y=1e-2, ls="--", c="gray")

ax[2].set_yscale("log")

# labels
ax[1].set_ylabel(L"$\mathcal{D}(\omega)$")
ax[2].set_xlabel(L"$\omega$")
ax[2].set_ylabel("error")

ax[2].legend(ncol=2)
plt.tight_layout()
Nz = length(params.z)
Λ = params.Λ
title = L"Mapping to a Wilson chain: $N_z = %$(Nz)$, $\Lambda = %$(Λ)$"
fig.suptitle(y=1.01, title)
display(fig)

##
using JLD2
