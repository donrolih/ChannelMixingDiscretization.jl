struct DiscretizationParams
    # NRG discretization parameter
    Λ::Float64
    # twisting parameters (always a vector, for one z input a lenght one vector)
    z::Vector{Float64}
    # Jsites
    J::Int64
    # number of samples for integration over ρ(ε(x))
    Nx::Int64
    # threshold to set an occasional small eigenvalue of ρ to zero
    threshold::Float64
    # discretization type
    disctype::Symbol
end

DiscretizationParams() = DiscretizationParams(2., 
                                              collect(range(0.05, 0.95, 10)),
                                              35,
                                              500000,
                                              1e-5,
                                              :adaptive
)

DiscretizationParams(Λ, zs, J) = DiscretizationParams(Λ, 
                                                      zs,
                                                      J,
                                                      500000,
                                                      1e-5,
                                                      :adaptive
)

DiscretizationParams(Λ, zs, J, disctype) = DiscretizationParams(Λ, 
                                                      zs,
                                                      J,
                                                      500000,
                                                      1e-5,
                                                      disctype
)

function logmesh(min, max, ratio)
    T = typeof(min)
    ωs = Array{T}(undef, 0)
    let ω = max
        while ω > min
            pushfirst!(ωs, ω)
            ω /= ratio
        end
    end
    return vcat(-reverse(ωs), ωs)
end

function logmesh(min, max, ratio, accpoint)
    T = typeof(min)
    ωs = Array{T}(undef, 0)
    let z = max
        rescale_factor = (max - accpoint)/max
        while z > min
            ω = accpoint + z*rescale_factor
            pushfirst!(ωs, ω)
            z /= ratio
        end
    end
    return vcat(-reverse(ωs), ωs)
end

function logmesh(mesh_min::BigFloat, mesh_max::BigFloat, mesh_ratio::BigFloat, mesh_accumulation::BigFloat, minatzero::BigFloat)
    posfreqs = BigFloat[]
    
    let ω = mesh_max
        while ω > mesh_min
            pushfirst!(posfreqs, ω)
            ω /= mesh_ratio
        end
    end
    
    posfreqs = vcat([minatzero], posfreqs)

    return vcat(-reverse(posfreqs), posfreqs)
end

# function logmesh(mesh_min, mesh_max, mesh_ratio, minatzero)
#     posfreqs = BigFloat[]
    
#     let ω = mesh_max
#         while ω > mesh_min
#             pushfirst!(posfreqs, ω)
#             ω /= mesh_ratio
#         end
#     end
    
#     posfreqs = vcat([minatzero], posfreqs)

#     return vcat(-reverse(posfreqs), posfreqs)
# end