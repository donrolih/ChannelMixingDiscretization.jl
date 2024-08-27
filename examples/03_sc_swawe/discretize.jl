using ChannelMixingDiscretization

# define the SC gap function

gap = 0.02
κ = 0.
Γ = 0.5

model = sWaveSC(gap, κ, Γ)

mesh_min = 1e-10
mesh_max = 1.
mesh_ratio = 1.01
mesh_accumulation = gap

# the input to the discretization function is a tabulated hybridization function
# frequencies
ωs = logmesh(mesh_min, mesh_max, mesh_ratio, mesh_accumulation)

# evaluate hybridization
ρs = hybridization.(ωs, Ref(model));

# number of points in the star representation (for positive AND negative branch)
# the length of the generated Wilson chain will be 2J
J = 30

# number of twisting parameters; equidistantly spaced from 1/Nz to 1
Nz = 4
zs = range(1/Nz, 1, Nz)

# grid type (see CPC paper)
gridtype = :fixed

# discretization parameter Λ for NRG
Λ = 2.
D = 1.

starH, chains = discretize(ωs, ρs,
                           mesh_min, mesh_max, mesh_ratio, mesh_accumulation,
                           J, zs, gridtype, Λ, D;
                           savechain=true, nrg_generatefolders=false)