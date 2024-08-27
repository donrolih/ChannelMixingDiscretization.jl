using ChannelMixingDiscretization

# hybridization function to be discretized
# Semicircular band with half-bandwidth D = 2t (Bethe lattice) 
function hybridization(ω; t=0.5)
    if abs(ω) < 2t
        return (1/(π*2t^2))*sqrt(4t^2 - ω^2)
    else
        return 0.
    end
end

mesh_min = 1e-20
mesh_max = 1.
mesh_ratio = 1.01

# no gap in this case
mesh_accumulation = 0.

# the input to the discretization function is a tabulated hybridization function
# frequencies
ωs = logmesh(mesh_min, mesh_max, mesh_ratio, mesh_accumulation)

# evaluate hybridization
ρs = hybridization.(ωs; t=0.5)

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