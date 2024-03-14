using DataInterpolations
using DelimitedFiles
using QuadGK
using Roots
includet("nrgutils.jl")

mutable struct DMFTSolver
    paramfile::String
    paramdict::Dict{String, String}
    # iteration number 
    iter::Int

    U::Float64
    occupancy_goal::Float64
    T::Float64

    # parameters of the mesh
    # in order to ensure stability of the discretization
    # they need to be of type BigFloat!
    mesh_max::BigFloat
    mesh_min::BigFloat
    mesh_ratio::BigFloat
    
    # immutable parameters
    # convergence criteria
    min_iter::Integer
    max_iter::Integer
    eps_prev::Float64
    eps_loc_imp::Float64
    eps_occuppancy::Float64
    
    solution_file::String
    checkpoint_file::String

    Delta_min::Float64

    # mixing parameter
    α ::Float64
    # number of cycles for adjusting the chemical potential
    max_mu_adjust::Integer
    mixing_method::Symbol

    ωs # ω mesh
    μ # current chemical potential
    Σ # current self-energy
    Δ # current hybridization function
    Gself # self-energy trick improved spectral function
    Gloc # local lattice Green's function
    Gloc_prev # previous iteration local lattice Green's function
    
    # incomplete inner constructor!
    function DMFTSolver(paramfile::String)
        # read the param file as a dictionary
        params = getparams(paramfile)
        
        U = parse(Float64, params["U1"])
        occupancy_goal = parse(Float64, params["occupancy_goal"])
        T = parse(Float64, params["T"])

        dos = params["nonintdos"]

        mesh_max = parse(BigFloat, params["broaden_max"])
        mesh_min = parse(BigFloat, params["broaden_min"])
        mesh_ratio = parse(BigFloat, params["broaden_ratio"])

        println(params["min_iter"])
        min_iter = parse(Int, params["min_iter"])
        max_iter = parse(Int, params["max_iter"])
        eps_prev = parse(Float64, params["eps_prev"])
        eps_loc_imp = parse(Float64, params["eps_loc_imp"])
        eps_occupancy = parse(Float64, params["eps_occupancy"])

        solution_file = params["solution_file"]
        checkpoint_file = params["checkpoint_file"]

        Delta_min = parse(Float64, params["Delta_min"])

        α = parse(Float64, params["mixing_alpha"])
        max_mu_adjust = parse(Int, params["max_mu_adjust"])

        ωs = convert.(Float64, logmesh(mesh_min, mesh_max, mesh_ratio))

        # initial Σ and μ
        if isfile(solution_file)
            # continue DMFT after failed convergence
            # TODO
        elseif isfile(checkpoint_file)
            # continue DMFT after interruption
            # TODO
        else
            # starting from scratch
            # we start with the Hartree shift
            μ = U/2.0
            Σ = zeros(length(ωs)) .+ U*occupancy_goal/2.0
        end
        return new(paramfile, params,
                   0, U, occupancy_goal, T,
                   mesh_max, mesh_min, mesh_ratio,
                   min_iter, max_iter, eps_prev, eps_loc_imp, eps_occupancy,
                   solution_file, checkpoint_file,
                   Delta_min, α, max_mu_adjust, :linear,
                   ωs, μ, Σ
        )
    end
end

function selfconsistency!(solver::DMFTSolver)
    ωs, μ, Σ = solver.ωs, solver.μ, solver.Σ
    solver.Gloc = ht.(ωs .+ μ .- Σ)
    Glocinv = inv.(solver.Gloc)
    solver.Δ = ωs .+ μ .- Σ .- Glocinv
    return nothing
end

function adjust_μ!(solver::DMFTSolver)
    old_μ = deepcopy(solver.μ)
    max_mu_adjust = solver.max_mu_adjust
    for i in 1:max_mu_adjust
        println("Adjusting μ at i = $(i)")
        update_μ!(solver::DMFTSolver)
        selfconsistency!(solver::DMFTSolver)
    end
    new_μ = solver.μ
    @info "Adjusted μ from $(old_μ) to $(new_μ)!"
    return nothing
end

"""
    Calculate occupancy for a given solver struct.
"""
function calc_occupancy(solver::DMFTSolver)
    ωs = solver.ωs
    Gtrial = calc_G(solver) # array 
    f = CubicSpline(-(1/pi) .* imag.(Gtrial), ωs) # be careful of the -(1/π) factor
    mesh_max = convert(Float64, solver.mesh_max) 
    occupancy, _ = quadgk(ω -> 2*f(ω)*expit(-ω/solver.T), -mesh_max, mesh_max)
    return occupancy
end

"""
    Calculate occupancy for a given hybridization, self-energy and chemical potential.
"""
function calc_occupancy(μ, Δ, Σ, ωs, T, mesh_max)
    Gtrial = 1 ./ (ωs .+ μ .- Δ - Σ)
    f = CubicSpline(-(1/pi) .* imag.(Gtrial), ωs) # be careful of the -(1/π) factor 
    occupancy, _ = quadgk(ω -> 2*f(ω)*expit(-ω/T), -mesh_max, mesh_max)
    return occupancy
end

function update_μ!(solver::DMFTSolver)
    μ, occupancy_goal, T = solver.μ, solver.occupancy_goal, solver.T
    println(μ)
    Δ, Σ, ωs = solver.Δ, solver.Σ, solver.ωs
    mesh_max = convert(Float64, solver.mesh_max)
    f(x) = calc_occupancy(x, Δ, Σ, ωs, T, mesh_max) - occupancy_goal
    sol = find_zero(f, μ)
    solver.μ = sol
    return nothing
end

"""
    Parse the parameter file. Returns a dictionary of param => value pairs.
"""
function getparams(file)
    # basically rewrite the parse_ini() function from ConfParser.jl
    # but modify it for my needs

    dict = Dict{String, String}()
    regex = r"^\s*([^=]*[^\s])\s*=\s*(.*)\s*$"
    for line in eachline(file)
        # skip comments and newlines
        occursin(r"^\s*(\n|\#|;)", line) && continue
        occursin(r"\w", line) || continue

        line = chomp(line)

        # parse key and value
        m = match(regex, line)
        if m !== nothing
            key::String, value::String = m.captures
            dict[key] = value
        end
    end
    return dict
end

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

"""
    Generate the Bethe lattice density of states on a logarithmic mesh.
"""
function bethe(min, max, ratio; t = big"0.5", savetofile=true, filename="Delta.dat", minvalue=big"1e-5")
    posfreqs = BigFloat[]
    values = BigFloat[]
    
    let ω = max
        while ω > min
            pushfirst!(posfreqs, ω)
            value = abs(ω) < 1. ? sqrt(4t^2 - ω^2)/(2pi*t^2) : minvalue
            pushfirst!(values, value)
            ω /= ratio
        end
    end

    # because the Bethe lattice is symmetric, we can just mirror for the negative values

    xs = vcat(-reverse(posfreqs), posfreqs)
    ys = vcat(reverse(values), values)

    if savetofile
        open(filename, "w") do io
            writedlm(io, [convert.(Float64, xs) convert.(Float64, ys)])
        end
    end

    return xs, ys
end

"""
    Calculate Green's function for a given solver object.
"""
function calc_G(solver::DMFTSolver)
    ωs = solver.ωs
    μ = solver.μ
    Δ = solver.Δ
    Σ = solver.Σ
    G = 1 ./ (ωs .+ μ .- Δ - Σ)
    return G
end

expit(x::Real) = inv(exp(-x) + one(x))

"""
    Hilbert transform of the Bethe lattice DOS.
"""
function ht(z; eps=1e-16)
    fixedim = imag(z) ≥ eps ? imag(z) : eps
    z = real(z) + im*fixedim
    return 2*(z - im*sign(imag(z))*sqrt(1 - z^2))
end

"""
    Given a path to the tabulated DOS returns an interpolation function.
"""
function interpolateDOS(path)
    if isfile(path)
        data = readdlm(path)
        return LinearInterpolation(data[:, 2], data[:, 1])
    else
        error("Could not find tabulated DOS!")
    end
end

function gfdiff(G₁, G₂, ωs, mesh_max)
    f₁ = CubicSpline(-(1/pi) .* imag.(G₁), ωs)
    f₂ = CubicSpline(-(1/pi) .* imag.(G₂), ωs)
    diff, _ = quadgk(ω -> (f₁(ω) - f₂(ω))^2, -mesh_max, mesh_max)
    return diff
end

# function fixhybridization!(solver::DMFTSolver)
#     Δ = solver.Δ
#     r = real.(Δ)
#     i = imag.(Δ)
#     Δ = r .+ im * ()
# end

struct Converged <: Exception end

struct FailedToConverge <: Exception end

"""
    Performs a step in the DMFT algorithm using an external program NRG Ljubljana. It does not return anything, but the output is written to a dat file. 
"""
function dmftstep!(solver::DMFTSolver)
    solver.iter += 1
    @info "Iteration $(solver.iter), min_iter $(solver.min_iter), max_iter $(solver.max_iter)"
    # fix hybridization function
    # fixhybridization!(solver::DMFTSolver)

    # now we have to solve the impurity problem
    # this should all happen in a subdirectory

    path = "imp_problem_iter=$(solver.iter)" # NRG is run here
    # NRG solver changes the self-energy
    nrg!(solver, path)
    
    # calculate impurity GF (self-energy trick improved)
    solver.Gself = calc_G(solver)
    
    # apply the self-consistency equation
    selfconsistency!(solver) # this sets new Δ and Gloc
    
    # Calculate the differences
    @info "Calculating differences and occupancies!"
    # difference between impurity Gself and local lattice GF Gloc
    mesh_max = convert(Float64, solver.mesh_max) # we do not need big floats for this
    @time diff_loc_imp = gfdiff(solver.Gself, solver.Gloc, solver.ωs, mesh_max)
    
    # difference between two consecutive local lattice GFs Gloc and Gloc_prev
    @time diff_prev = gfdiff(solver.Gloc, solver.Gloc_prev, solver.ωs, mesh_max)

    solver.Gloc_prev = deepcopy(solver.Gloc)  # deepcopy is maybe a bit paranoic

    # calculate occupancy
    @time occupancy = calc_occupancy(solver)
    
    # calculate the difference between occupancy and the target occupancy
    diff_occupancy = abs(occupancy - solver.occupancy_goal)

    # maybe join all of these diff into a struct "DMFTMonitor"
    # which can be saved at each iteration
    # check for convergence
    eps_loc_imp = solver.eps_loc_imp
    eps_prev = solver.eps_prev
    eps_occupancy = solver.eps_occuppancy
    min_iter = solver.min_iter
    max_iter = solver.max_iter

    @info "Occupancy: $(occupancy)" 
    @info "Difference between the local GF and impurity GF: $(diff_loc_imp)"
    @info "Difference between two consecutive local GFs: $(diff_prev)"
    @info "Difference between occupancy and target occupancy: $(diff_occupancy)"

    if diff_loc_imp < eps_loc_imp && diff_prev < eps_prev && diff_occupancy < eps_occupancy && solver.iter ≥ min_iter
        # we have converged
        # save
        mkpath("solution/")
        open("solution/Gloc.dat", "w") do io
          writedlm(io, [real.(solver.ωs) real(solver.Gloc) imag.(solver.Gloc)])
        end
        ##
        open("solution/Sigma.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Σ) imag.(solver.Σ)])
        end

        open("solution/Gself.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Gself) imag.(solver.Gself)])
        end

        open("solution/Delta.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Δ) imag.(solver.Δ)])
        end
        # throw the error
        throw(Converged)
    elseif solver.iter == max_iter
        mkpath("solution/")
        
        open("solution/Gloc.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real(solver.Gloc) imag.(solver.Gloc)])
        end

        open("solution/Sigma.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Σ) imag.(solver.Σ)])
        end

        open("solution/Gself.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Gself) imag.(solver.Gself)])
        end

        open("solution/Delta.dat", "w") do io
            writedlm(io, [real.(solver.ωs) real.(solver.Δ) imag.(solver.Δ)])
        end
        
        throw(FailedToConverge)
    end
    
    # TODO: adjust occupancy
    @info "Adjusting μ!"
    adjust_μ!(solver::DMFTSolver)

    @info "End of DMFT step at iteration $(solver.iter)!"
end

function solve!(solver::DMFTSolver)
    mixing_method = solver.mixing_method
    if mixing_method == :linear
        while true # we exit by triggering an Exception in dmftstep! function
            Δ_in = deepcopy(solver.Δ)
            dmftstep!(solver)
            solver.Δ = solver.α .* solver.Δ + (1 - solver.α) .* Δ_in
        end
    else
        # TODO: use the Broyden method from NonlinearSolve.jl
        error("currently only linear mixing is implemented!")
    end
    nothing
end

# This is how we run (catching the exceptions)
# try
#     solve!(solver)
# catch e
#     if isa(e, Converged)
#         @info "We have converged!"
#     elseif isa(e, FailedToConverge)
#         @warn "We have reached maximum iteration number without converging!"
#     end
# end