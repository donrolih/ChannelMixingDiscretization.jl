using ChannelMixingDiscretization
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

# saving
savechains(chains);