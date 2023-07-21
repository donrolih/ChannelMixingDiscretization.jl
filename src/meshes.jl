struct DiscretizationParams
    # NRG discretization parameter
    Λ::Real
    # twisting parameters (always a vector, for one z input a lenght one vector)
    z::Vector
    # Jsites
    J::Integer
    # number of samples for integration over ρ(ε(x))
    Nx::Integer
    # threshold to set an occasional small eigenvalue of ρ to zero
    threshold::Real
end

DiscretizationParams() = DiscretizationParams(2., 
                                              collect(range(0.05, 0.95, 10)),
                                              35,
                                              500000,
                                              1e-5
)

DiscretizationParams(Λ, zs, J) = DiscretizationParams(Λ, 
                                                      zs,
                                                      J,
                                                      500000,
                                                      1e-5
)

"""
    Mesh super type.
"""
abstract type Mesh end

"""
    Logarithmic mesh.
"""
struct LogMesh <: Mesh
    ω0::Real
    Nω::Integer
    D::Real
end

LogMesh() = LogMesh(1e-12, 5000, 1.)

"""
    Linear mesh.
"""
struct LinMesh <: Mesh
    ω0::Real
    Nω::Integer
    D::Real
end