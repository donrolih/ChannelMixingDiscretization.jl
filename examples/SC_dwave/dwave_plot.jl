##
using ChannelMixingDiscretization
##
mesh = LogMesh()

# energy gap
Δ = 0.1
# particle-hole symmetry
κ = 1.

model = sWaveSC(Δ, 0., κ, 0., 0.5/pi)

Nz = 64
J = 35

zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
params = DiscretizationParams(2., zs, J, 500000, 1e-5);

# hybridization function
function ρ(ω; Δ=0.1, t0=1., t2=0.1)
    σ0 = [1. 0.; 0. 1.]
    σ1 = [0. 1.; 1. 0.]
    
    d0 = (t0^2 + 2t2^2)*abs(ω)/Δ
    d1 = t0*t2*ω*abs(ω)/Δ^2
    return (d0 .* σ0) .+ (d1 .* σ1) 
end

##
ωs = generateω(mesh, model)
starH = discmodel(ρ, mesh, model, params);
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
pos = ωs .> 0
for i in 1:2
    ax[1].plot(ωs[pos], pauli_coeffs[i, pos], label=labels[i], color="C$(i-1)")
    ax[1].plot(ωs[pos], pauli_coeffs_re[i, pos], ".", label=labels[i], color="C$(i-1)")
    difpos = abs.(pauli_coeffs_re[i, pos] .- pauli_coeffs[i, pos])
    ax[2].plot(ωs[pos], difpos, label=labels[i], color="C$(i-1)")

    neg = .!pos
    ax[1].plot(ωs[neg], pauli_coeffs[i, neg], color="C$(i-1)")
    ax[1].plot(ωs[neg], pauli_coeffs_re[i, neg], ".", color="C$(i-1)")
    difneg = abs.(pauli_coeffs_re[i, neg] .- pauli_coeffs[i, neg])
    ax[2].plot(ωs[neg], difneg, color="C$(i-1)")
end

# ax[1].set_ylim(-0.2, 0.2)

ax[2].axhline(y=1e-2, ls="--", c="gray")
ax[2].axvline(x=-Δ, ls="--", c="gray")
ax[2].axvline(x=Δ, ls="--", c="gray")

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
pos = ωs .> 0
for i in 1:2
    ax[1].plot(ωs[pos], pauli_coeffs[i, pos], label=labels[i], color="C$(i-1)")
    ax[1].plot(ωs[pos], pauli_coeffs_chain[i, pos], ".", label=labels[i], color="C$(i-1)")
    difpos = abs.(pauli_coeffs_chain[i, pos] .- pauli_coeffs[i, pos])
    ax[2].plot(ωs[pos], difpos, label=labels[i], color="C$(i-1)")

    neg = .!pos
    ax[1].plot(ωs[neg], pauli_coeffs[i, neg], color="C$(i-1)")
    ax[1].plot(ωs[neg], pauli_coeffs_chain[i, neg], ".", color="C$(i-1)")
    difneg = abs.(pauli_coeffs_chain[i, neg] .- pauli_coeffs[i, neg])
    ax[2].plot(ωs[neg], difneg, color="C$(i-1)")
end

# ax[1].set_ylim(-0.2, 0.2)

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
