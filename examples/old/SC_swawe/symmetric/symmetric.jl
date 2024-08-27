##
using ChannelMixingDiscretization
##
mesh = LogMesh()

# energy gap
Δ = 0.1
# particle-hole asymmetry
κ = 0.

model = sWaveSC(Δ, 0., κ, 0., 0.5/pi)

Nz = 32
J = 35

zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
params = DiscretizationParams(2., zs, J, 500000, 1e-5);
##
ωs = generateω(mesh, model)
ρ = hybri(model, mesh);
##
# mapping to star model
starH = discmodel(ρ, mesh, model, params);
# mapping to chain
chains = maptochains(starH; m=J);
##
# saving
savechains(chains)