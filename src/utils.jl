########################
#GENERATING DATA
########################

"""
    Generate a vector of Nz twisting parameters, linearly distributed on (0, 1).
"""
function twistingparameters(Nz::Integer)
    return collect(range(0.5/Nz, 1 - 0.5/Nz, Nz))
end

"""
    Generates ω support of ρ for a given mesh type and model. 
"""
function generateω(mesh::LogMesh, model::PhysicalModel; offset=1e-30)
    # mesh parameters
    ω0, Nω, D = mesh.ω0, mesh.Nω, mesh.D
    Δ = model.Δ

    negfreqs = exp.(range(log(ω0), log(D - Δ), Int(Nω/2) - 1)) .+ Δ
    posfreqs = exp.(range(log(ω0),log(D - Δ), Int(Nω/2) - 1)) .+ Δ

    return vcat(-reverse(negfreqs), [-Δ - offset, Δ + offset], posfreqs)
end

########################
# DATA ANALYSIS
########################
"""
    Coefficients of the expansion of the D ∈ su(2) in the basis of Pauli matrices.
"""
function paulibasis(D::Matrix)
    σ1 = [0. 1.; 1. 0.]
    σ2 = [0. -im; im 0.]
    σ3 = [1. 0; 0. -1.]

    coeffs = zeros(ComplexF64, 4)
    coeffs[1] = tr(D)
    coeffs[2] = tr(σ1*D)
    coeffs[3] = tr(σ2*D)
    coeffs[4] = tr(σ3*D)

    return coeffs./2
end

"""
    Get the Pauli coefficients of a hybridization function for all ω.
"""
function getpaulicoeffs(hybri)
    Nω = size(hybri, 1)
    allcoeffs = zeros(4, Nω)
    for i in 1:Nω
        D = hybri[i, :, :]
        coeffs = paulibasis(D)
        allcoeffs[:, i] = real.(coeffs)
    end
    return allcoeffs
end