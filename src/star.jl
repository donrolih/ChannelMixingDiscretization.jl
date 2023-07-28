struct StarHamiltonian 
    T::Dict
    E::Dict
    zs::Vector
end

struct Ticker
    R::AbstractArray 
    iR::Interpolations.Extrapolation
    function Ticker(ωs::AbstractArray, weights::AbstractArray)
        R = cumul_integrate(ωs, weights)
        iR = linear_interpolation(R, ωs, extrapolation_bc=Line())
        return new(R, iR)
    end
end

# guiding function, also called ν(x)
function ε(x::Real,
           ticker::Ticker,
           discparams::DiscretizationParams,
           mesh::Mesh)
    Λ = discparams.Λ
    R = ticker.R
    RD = last(R)
    iR = ticker.iR
    if x <= 2
        return mesh.D
    else
        # adaptive
        return iR(RD*Λ^(2 - x))
        # fixed log-gap
        # return Λ^(2 - x)
    end
end

function ε(xs::AbstractArray,
           ticker::Ticker,
           discparams::DiscretizationParams,
           mesh::Mesh)
    return [ε(x, ticker, discparams, mesh) for x in xs]
end

function getweights(ρs::Vector)
    el = first(ρs)
    if isa(el, Number) || isa(el, Complex)
        return ρs
    elseif isa(el, Matrix)
        # compute the squared Hilbert-Schmidt norm (it probably should not be squared, but I am following their implementation)
        # these are the weights for the adaptive scale
        println("Getting weights (matrices) ...")
        norms = [norm(ρ) for ρ in ρs]
        return norms
    else
        error("something is wrong with the dimensions of the hybridisation function!")
    end
end

"""
    Returns the T(x) and E(x) functions with a given ρ.
"""
function getTEfunctions(ωs::Vector,
                        ρs::Vector,
                        params::DiscretizationParams,
                        mesh::Mesh,
                        weights::Vector)

    ωbranches = Dict()
    ρbranches = Dict()
    weightbranches = Dict()

    positivemask = ωs .> 0
    ωbranches[1] = ωs[positivemask]
    ωbranches[-1] = reverse(ωs[.!positivemask])

    ρbranches[1] = ρs[positivemask]
    ρbranches[-1] = reverse(ρs[.!positivemask])

    # weights = getweights(ρs)
    weightbranches[1] = weights[positivemask]
    weightbranches[-1] = reverse(weights[.!positivemask])

    Tfunctions = Dict()
    Efunctions = Dict()
    signs = [-1, 1]

    for sign in signs
        # select the positive/negative branch
        ρbranch = ρbranches[sign]
        ωbranch = ωbranches[sign]
        
        weightbranch = weightbranches[sign]
        ticker = Ticker(sign*ωbranch, weightbranch)

        # get the R(ω) function, see Appendix 1 of PRB 93, 035102 (2016)
        integratedρ = sign .* cumul_integrate(ωbranch, ρbranch)
        Rfunc = linear_interpolation(sign.*ωbranch, integratedρ, extrapolation_bc=Line())
        
        function Tfunc(x)
            lb = ε(x + 1, ticker, params, mesh)
            ub = ε(x, ticker, params, mesh)
            arg = Rfunc(ub) - Rfunc(lb)
            if arg >= 0.
                return sqrt(arg)
            elseif abs(arg) < 1e-6
                return 0.
            else
                println(arg)
                error("negative argument when calculating the hopping function!")
            end
        end


        Tfunctions[sign] = Tfunc
        # it has to be J + 2 because of the extended bound of integration below
        xs = range(1, params.J + 2, params.Nx)
        # inverse of the R function
        iRfunc = linear_interpolation(integratedρ, sign.*ωbranch, extrapolation_bc=Line())

        # change the integration order to avoid big+small
        RR = vcat(cumul_integrate(reverse(xs), Rfunc.(ε(reverse(xs), ticker, params, mesh))))
        RRfunc = linear_interpolation(xs, reverse(RR), extrapolation_bc=Line())
        Efunc(x) = sign*(iRfunc.(RRfunc(x+1) - RRfunc(x)))
        Efunctions[sign] = Efunc
    end

    return Tfunctions, Efunctions
end

function getTElists(Efunctions, Tfunctions, ρ, params, Nbands)
    zs = params.z
    J = params.J
    
    Ts = Dict()
    Es = Dict()
    if Nbands == 1
        for sign in [-1, 1]
            E = zeros(Complex{BigFloat}, J, length(zs))
            T = zeros(Complex{BigFloat}, J, length(zs))
            for j in 1:J
                for (k, z) in enumerate(zs)
                    x = big(j + z)
                    ϵ = Efunctions[sign](x)
                    t = Tfunctions[sign](x)
                    E[j, k] = ϵ
                    T[j, k] = t
                end
            end
            Es[sign] = E
            Ts[sign] = T
        end
    else
        for sign in [-1, 1]
            # shape of E and T is: (number of J sites, Nz, Nbands, Nbands)
            E = zeros(Complex{BigFloat}, J, length(zs), Nbands, Nbands)
            T = zeros(Complex{BigFloat}, J, length(zs), Nbands, Nbands)
            for j in 1:J
                for (k, z) in enumerate(zs)
                    x = big(j + z)
                    ϵ = [Efunctions[sign][i](x) for i in 1:Nbands]
                    ϵ = diagm(ϵ)
                    E[j, k, :, :] = ϵ
                    ϵiold = Inf
                    Ux = zeros(Complex{BigFloat}, Nbands)
                    for i in 1:Nbands
                        ϵi = ϵ[i, i]
                        if abs(ϵiold - ϵi) > 1e-12
                            M = Hermitian(ρ(ϵi))
                            Ui = eigvecs(M)[:, i]
                            ϵiold = ϵi
                        else
                            println("Found a degeneracy!")
                            if i == 1 error("impossible!") end
                            # we have to choose the same vector as before
                            Ui = Ux[:, end]
                        end
                        Ux = hcat(Ux, Ui)
                    end
                    Ux = Ux[:, 2:end]
                    tfuncs = Tfunctions[sign]
                    Td = zeros(Complex{BigFloat}, Nbands)
                    for i in 1:Nbands
                        Tdi = tfuncs[i](x).*Ux[:, i]
                        Td = hcat(Td, Tdi)
                    end
                    T[j, k, :, :] = Td[:, 2:end]'
                end
            end
            Es[sign] = E
            Ts[sign] = T
        end
    end

    return Ts, Es
end

function discmodel(ρ::Function, mesh::Mesh, model::PhysicalModel, params::DiscretizationParams)
    ωs = generateω(mesh, model)
    ρs = [ρ(ω) for ω in ωs]
    Tfunctions = Dict()
    Efunctions = Dict()
    weights = getweights(ρs)
    if isa(first(ρs), Number)
        # single band case
        Nbands = 1
        println("Single band case!")
        println("Getting functions of representative energies and hoppings ...")
        Tfunctions, Efunctions = getTEfunctions(ωs, ρs, params, mesh, weights)
        println("Done!")
        println("Generating discrete model ...")
        Ts, Es = getTElists(Efunctions, Tfunctions, ρ, params, Nbands)
        println("Done!")
        return StarHamiltonian(Ts, Es, params.z)
    elseif isa(first(ρs), Matrix)
        # multiband
        Nbands, _ = size(first(ρs))
        println("Multi band case: Nbands = $(Nbands)!")
        println("Diagonalizing the hybridisation function ...")
        ρeval = [eigvals(M) for M in ρs]
        ρeval = reduce(vcat, transpose.(ρeval))
        println("Done!")
        # _, Nbands = size(ρeval)
        # set an occasional small eigenvalue to zero

        println("Getting functions of representative energies and hoppings ...")
        Tplus = Vector{Function}(undef, Nbands)
        Tmin = Vector{Function}(undef, Nbands)
        Eplus = Vector{Function}(undef, Nbands)
        Emin = Vector{Function}(undef, Nbands)
        for i in 1:Nbands
            println("Calculating for band number $(i) ...")
            start = time()
            tfs, efs = getTEfunctions(ωs, ρeval[:, i], params, mesh, weights)
            stop = time()
            elapsed = stop - start
            println("Time elapsed for band number $(i): $(round(elapsed; digits=3)) s")
            Eplus[i] = efs[1]
            Emin[i] = efs[-1]
            Tplus[i] = tfs[1]
            Tmin[i] = tfs[-1]
        end
        Efunctions[1] = Eplus
        Efunctions[-1] = Emin
        Tfunctions[1] = Tplus
        Tfunctions[-1] = Tmin
        println("Done!")

        println("Generating discrete model ...")
        Ts, Es = getTElists(Efunctions, Tfunctions, ρ, params, Nbands)
        println("Done!")
        return StarHamiltonian(Ts, Es, params.z)
    else
        error("something is wrong with the dimensions of the hybridisation function!")
    end
end