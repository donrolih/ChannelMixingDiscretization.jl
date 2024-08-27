using ChannelMixingDiscretization
##
# offset is important (ω0!!! if it is too small, tou get divergences at zero!)
mesh = LogMesh(1e-4, 1000, 1.)
model = sWaveSC(0., 0., 0., 0., 0.5/pi)
Nz = 64
zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
J = 35
params = DiscretizationParams(2., zs, J, 500000, 1e-5);
##
function ρ(ω::Real)
    σ0 = [1. 0; 0. 1.]
    σ1 = [0. 1.; 1. 0.]
    σ2 = [0. -im ; im 0.]
    return σ0 + 0.4*σ1 + 0.5*(0.1 + ω^2)*σ2
end;

##
# mapping to star model
starH = discmodel(ρ, mesh, model, params);
# mapping to chain
chains = maptochains(starH; m=J);
##
# saving
savechains(chains)