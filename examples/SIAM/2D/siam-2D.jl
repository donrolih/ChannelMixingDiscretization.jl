##
using ChannelMixingDiscretization
##
mesh = LogMesh(1e-4, 2500, 1.)
model = sWaveSC(0., 0., 0., 0., 0.5/pi)

Nz = 32
zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
J = 35

params = DiscretizationParams(2., zs, J, 500000, 1e-5);