using ChannelMixingDiscretization
using DelimitedFiles
##
# offset is important (ω0!!! if it is too small, you get divergences at zero!)
mesh = LogMesh()
model = Flat()

Nz = 32
J = 60
# zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
zs = collect(range(1, Nz, Nz)/Nz)
params = DiscretizationParams(1.8, zs, J, 500000, 1e-5);

# hybridization function to map
ρ = hybri(model, mesh);
##
# mapping to star model
starH = discmodel(ρ, mesh, model, params);
##
# mapping to Wilson chain
chains = maptochains(starH; m=J);
##
# spectral function
freqs, weights = spectralfunction(chains);

##
# saving
# savechains(chains)