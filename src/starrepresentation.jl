# NOTE:
# This file implements the discretization of a given hybridization function
# The end result is given as a struct StarHamiltonian, which holds
# hopping matrices T and on-site energy matrices E for each twisting parameter z

# N. B.: we use BigFloats throughout; this improves acccuracy when hybridization
# has divergences

# Some references:
#   - R. Zitko, CPC 180 (2009), 1271-1276,
#   - Liu et.al., PRB 93, 035102 (2016).

################################################################################
# GRID STRUCTS and GUIDING FUNCTIONS FOR FIXED AND ADAPTIVE GRIDS
################################################################################

struct AdaptiveGrid
    interp_iR::DataInterpolations.AbstractInterpolation
    Rtot::BigFloat
    function AdaptiveGrid(ωs, weights)
        R = cumintegrate(ωs, weights)
        Rtot = R[end]
        interp_iR = LinearInterpolation(ωs, R)
        return new(interp_iR, Rtot)
    end
end

struct FixedGrid
    Λ::BigFloat  # discretization parameter
    D::BigFloat  # half-bandwidth, the spectal fn. is contained in [-D, D]
    Δ::BigFloat  # gap (this changes the accumulation point of the guiding fn)
end

FixedGrid() = FixedGrid(big"2.", big"1.", big"0.0")

function makeguiding(grid::AdaptiveGrid)
    iR = grid.interp_iR
    Rtot = grid.Rtot
    return x -> x ≤ big"2." ? iR(Rtot) : iR(Rtot*(big"2." ^ (big"2." - x)))
end

function makeguiding(grid::FixedGrid)
    return x -> x ≤ big"2." ? grid.D : (grid.D - grid.Δ)*grid.Λ^(big"2." - x) + grid.Δ
end

################################################################################
# DISCRETIZER STRUCT
################################################################################

struct Discretizer
    ν               # guiding function
    interp_R        # interpolation of the ωs -> weights function
    interp_iR       # interpolation of the weights --> ωs function
    interp_Rint     # interpolation of x -> ∫R(ν(x'))dx'; integral 1 to x
    bsign           # sign of ωs (positive or negative branch)
    function Discretizer(ωbranch, ρbranch, weightbranch, 
                         bsign, gridtype, xmax, Nx)
        grid = gridtype == :fixed ? FixedGrid() : AdaptiveGrid(bsign .* ωbranch, weightbranch)
        ν = makeguiding(grid)
        R = cumintegrate(bsign .* ωbranch, ρbranch)
        interp_R = LinearInterpolation(R, bsign .* ωbranch)
        interp_iR = LinearInterpolation(bsign .* ωbranch, R)
        xs = range(one(BigFloat), xmax, Nx)
        Rint = cumintegrate(xs, interp_R(ν.(xs)))
        interp_Rint = LinearInterpolation(Rint, xs)
        return new(ν, interp_R, interp_iR, interp_Rint, bsign)
    end
end

Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype) = Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype, big"200.", 10000)

function w(disc::Discretizer, x::BigFloat)
    ν = disc.ν
    interp_R = disc.interp_R
    upperbound = ν(x)
    lowerbound = ν(x + one(BigFloat))
    return interp_R(upperbound) - interp_R(lowerbound)
end

function w(disc::Discretizer, xs::Vector{BigFloat})
    ws = Vector{BigFloat}(undef, size(xs))
    for (i, x) in enumerate(xs)
        ws[i] = w(disc, x) 
    end
    return ws
end

function ε(disc::Discretizer, x::BigFloat)
    interp_iR = disc.interp_iR
    interp_Rint = disc.interp_Rint
    bsign = disc.bsign
    return bsign .* interp_iR(interp_Rint(x + one(BigFloat)) - interp_Rint(x))
end

function ε(disc::Discretizer, xs::Vector{BigFloat})
    εs = Vector{BigFloat}(undef, size(xs))
    for (i, x) in enumerate(xs)
        εs[i] = ε(disc, x)
    end
    return εs
end

function maskrho(ρeval::Vector{BigFloat}, bsign, mask)
    ρbranch = bsign == 1 ? ρeval[mask] : reverse(ρeval[mask])
    return [ρbranch]
end

function maskrho(ρeval::Matrix{BigFloat}, bsign, mask)
    # number of bands is the number of ρeval columns
    N = size(ρeval, 2)
    ρbranches = Vector{Vector{BigFloat}}(under, N)
    for (i, rhovec) in enumerate(eachcol(ρeval))
        ρbranches[i] = bsign == 1 ? rhovec[mask] : reverse(rhovec[mask])
    end
    return ρbranches
end

function calculatediscretizers(ωs::Vector{BigFloat}, ρeval::VecOrMat{BigFloat}, weights::Vector{BigFloat}, bsign::Integer, gridtype::Union{FixedGrid, AdaptiveGrid})
    @assert (bsign == 1) | (bsign == -1) "the value of bsign is not 1 or -1"
    
    mask = bsign == 1 ? ωs .> zero(BigFloat) : ωs .< zero(BigFloat)

    ωbranch = bsign == 1 ? ωs[mask] : reverse(ωs[mask])
    weightbranch = bsign == 1 ? weights[mask] : reverse(weights[mask])

    # always returns a vector, even if ρeval has one band
    # length of ρbranches corresponds to number of bands
    # this is the only difference between different bands
    # CAREFUL: the input ρeval should already be of shape (length(ωs), Nbands)
    # Each column of this matrix represents the positive eigenvalues of the hybridization matrix
    # The eigenvalues are positive because this is a positive semidefinite matrix
    ρbranches = maskrho(ρeval, bsign, mask) 
    discretizers = Vector{Discretizer}(undef, length(ρbranches))
    for (i, ρbranch) in enumerate(ρbranches)
        discretizers[i] = Discretizer(ωbranch, ρbranch, weightbranch, bsign, gridtype)
    end
    return discretizers
end

################################################################################
# EVALUATE COEFFICIENTS for DIFFERENT DISCRETIZER STRUCTS
################################################################################
# N.B.: we use multiple dispatch to differentiate between one channel and multichannel cases

# one-channel version
function evaluatecoefficients(discpos::Discretizer, discneg::Discretizer, J, zs)
    # We store the matrices in a dictionary
    Es = Dict()
    Ts = Dict()
    
    # On-site energy matrices and hopping matrices (here, they are numbers)
    Es[-1], Ts[-1] = evaluatebranch(discneg, J, zs)
    Es[1], Ts[1] = evaluatebranch(discpos, J, zs)

    return Es, Ts
end

# multichannel version (requires the hybridization function)
function evaluatecoefficients(discpos::Vector{Discretizer}, discneg::Vector{Discretizer}, J, zs, ρ)
    # We store the matrices in a dictionary
    Es = Dict()
    Ts = Dict()

    Es[-1], Ts[-1] = evaluatebranch(discneg, J, zs, ρ)
    Es[1], Ts[1] = evaluatebranch(discpos, J, zs, ρ)

    return Es, Ts
end

# one-channel version
function evaluatebranch(disc::Discretizer, J, zs)
    # one channel
    # number of twisting parameters
    Nz = length(zs)
    # we are storing numbers for each x = j + z
    E = zeros(Complex{BigFloat}, J, Nz)
    T = zeros(Complex{BigFloat}, J, Nz)
    for j in 1:J
        for (k, z) in enumerate(zs)
            x = big(j + z)  # evaluate the representative en. at this point
            ϵ = ε(disc[1], x)  # representative energy
            t = sqrt(w(disc[1], x))
            E[j, k] = ϵ
            T[j, k] = t
        end
    end
    return E, T
end

# multi-channel version (requires the hybridization function)
function evaluatebranch(disc::Vector{Discretizer}, J, zs, ρ)
    # multichannel case
    Nbands = length(disc)
    # now we store matrices for each x = j + z
    E = zeros(Complex{BigFloat}, J, Nz, Nbands, Nbands)
    T = zeros(Complex{BigFloat}, J, Nz, Nbands, Nbands)
    for j in 1:J
        for (k, z) in enumerate(zs)
            x = big(j + z)
            # we evaluate now for each band
            ϵ = [ε(disc[i], x) for i in 1:Nbands]
            # On-site energy matrices are diagonal in the star rep.
            # Construct diagonal matrix
            ϵ = diagm(ϵ)
            E[j, k, :, :] = ϵ
            # Now comes the evaluation of hopping matrices
            # Here more care is required
            # We have to diagonalize the evaluation of ρ(ϵ[i, i])
            # ρ(ϵi) is Hermitian (even semi positive definite) since
            # it is the hybridization function
            # Q: can the semi-definitness be exploited? (currently it isn't)
            # Maybe we can exploit the fact that all eigenvalues are nonnegative?
            ϵiold = Inf
            # this will store the eigenvectors
            Ux = zeros(Complex{BigFloat}, Nbands)
            for i in 1:Nbands
                ϵi = ϵ[i, i]
                if abs(ϵiold - ϵi) > 1e-100
                    # eigenvalues are different (we use big floats)
                    M = Hermitian(ρ(ϵi))
                    Ui = eigvecs(M)[:, i]
                    ϵiold = ϵi
                else
                    @info "Found a degeneracy for z-number z = $(z) at index J = $(J)"
                    # we have to choose the same eigenvector as before
                    Ui = Ux[:, end]
                end
            end
            Ux = Ux[:, 2:end]
            Td = zeros(Complex{BigFloat}, Nbands)
            for i in 1:Nbands
                Tdi = sqrt(w(disc[i], x)) .* Ux[:, i]
                Td = hcat(Td, Tdi)
            end
            T[j, k, :, :] = Td[:, 2:end]'
        end
    end
    return E, T
end