function gramschmidt(u, Q; maxiter=3, α=big"0.5")
    # println(size(u'*u))
    r = norm(u'*u)
    r1 = big"0."
    u_start = u
    for i in 1:maxiter
        u = u - Q*Q'*u
        r1 = norm(u'*u)
        # println("pre norm: ", r)
        # println("after norm: ", r1)
        if r1 > α*r
            break
        end
        r = r1
    end
    if r1 <= α*r
        println("WARNING! Loss of orthogonality at Gram-Schmidt!")
    end

    return u
end

"""
    Given a symmetric matrix A and a starting (block) vector q tridiagonalize the matrix A where the basis consists of Krylov vectors constructed from q.

    Implementation follows:
        - Grimes et. al., A shifted block Lanczos algorithm for solving sparse symmetric generalized eigenproblems, SIAM J. Matrix Anal. Applied 15, 228-272, 1994.
"""
function tridiagonalize(A, q, m)
    # A has shape 2JI x 2JI, q is a 'thin matrix' of shape 2JI x I
    # I is (this is also Nbands)
    n = size(q, 2)
    # set the number of iterations as 2J
    m = m === nothing ? Int(size(A)[1]/n) : m
    Al = zeros(Complex{BigFloat}, m, n, n)
    B = zeros(Complex{BigFloat}, m+1, n, n)
    Bu = zeros(Complex{BigFloat}, m, n, n)
    Q = hcat(zeros(BigFloat, size(q)), q)
    # println("starting vector ortho: ", float.(q'*q))
    for j in 1:m
        Qj = Q[:, end-(n-1):end]
        Uj = A*Qj - Q[:, (end - 2n + 1):(end - n)]*(B[j, :, :]')
        Aj = Qj'*Uj
        Rj = Uj - Qj*Aj
        F = qr(Rj)
        # new Q
        newQj = Matrix(F.Q)

        # reorthogonalize
        newQj = gramschmidt(newQj, Q[:, n+1:end])
        Q = hcat(Q, newQj)
        
        if (j != m) && norm(F.R, 1) < 1e-20
            println("WARNING! Bad Krylov space!")
        end
        
        # store the results
        Al[j, :, :] = Aj
        B[j + 1, :, :] = Matrix(F.R)
        Bu[j, :, :] = Matrix(F.R)'
    end
    Bl = B[2:end-1, :, :]
    # return a dictionary of arrays; keys denote the (off)-diagonal index
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end

"""
    This is an implementation of the numericaly stable Lanczos algorithm called
    Rutihauser-Kahan-Pal-Walker algorithm.
    Using this algorithm one can use double precision for the mapping to a 
    tight-binding (Wilson) chain.

    See the papers: 
        - Gragg and Harrod, Numer. Math 44, 317-335 (1984),
        - Gautschi,  Journal of Computational and Applied Mathematics 178, 215–234 (2005),
        - de Vega et. al., Phys. Rev. B 92, 155126 (2015) (mentioned in the Appendix).

    Notation follows Gragg and Harrod's paper.
"""
function rkpw(N, nodes, weights)
    @assert N > 0 "N should be positive!"
    Nmax = size(nodes, 1)
    @assert N ≤ size(nodes, 1) "N cannot be greater than the number of nodes"
    
    p0 = deepcopy(nodes)
    p1 = zeros(Nmax)
    p1[1] = weights[1]

    for n in 1:(Nmax-1)
        pn = weights[n + 1]
        γ = 1.
        σ = 0.
        τ = 0
        λ = nodes[n+1]
        for k in 1:n+1
            ρ = p1[k] + pn
            tmp = γ*ρ
            τσ = σ
            if ρ ≤ 0
                γ = 1
                σ = 0
            else
                γ = p1[k]/ρ
                σ = pn/ρ
            end
            τk = σ*(p0[k] - λ) - γ*τ
            p0[k] = p0[k] - (τk - τ)
            τ = τk
            if σ ≤ 0 
                pn = τσ*p1[k]
            else
                pn = τ^2/σ
            end
            τσ = σ
            p1[k] = tmp
        end
    end
    return [p0[1:N] p1[1:N]]
end

function blocklanczos(H::Matrix{Complex{T}}, Q1::VecOrMat{Complex{T}}, J::Integer) where T <: Real
    @assert Q1'*Q1 ≈ I "the initial block vector for Lanczos tridiagonalization should be an orthogonal matrix"
    n, p = size(Q1)
    @assert J ≤ n/p "the number of steps is too large for matrix size" 
    # the storage of Q blocks
    Q = Vector{Matrix{T}}(undef, J)
    Q[1] = Q1
    A = Vector{Matrix{T}}(undef, J)
    B = Vector{Matrix{T}}(undef, J)
    B[1] = zeros(T, p, p)
    for j in 1:J
        Qjj = j == 1 ? zeros(T, n, p) : Q[j - 1]
        U = H*Q[j] - Qjj * B[j]'
        A[j] = Q[j]'* U
        R = U - Q[j]*A[j]
        F = qr(R)
        if j == J break end
        Q[j+1] = Matrix(F.Q)
        B[j+1] = Matrix(F.R)
    end
    # set the last diagonal element
    return A, B[2:end], adjoint.(B[2:end])
end

function tridiagonalizeONEBAND(A, q, m)
    n = size(q, 2)
    # set the number of iterations as 2J
    m = m === nothing ? Int(size(A)[1]/n) : m
    # Ss = zeros(Float64, m, n, n)
    Al = zeros(Complex{BigFloat}, m, n, n)
    B = zeros(Complex{BigFloat}, m+1, n, n)
    Bu = zeros(Complex{BigFloat}, m, n, n)
    Q = hcat(zeros(BigFloat, size(q)), q)
    for j in 1:m
        Qj = Q[:, end-(n-1):end]
        Uj = A*Qj - Q[:, (end - 2n + 1):(end - n)]*(B[j, :, :]')
        Aj = Qj'*Uj
        Rj = Uj - Qj*Aj
        # ortho = dot(Rj[:, 1], Rj[:, 2])
        # println("are the columns of Rj orthogonal? ", norm(ortho) < 1e-20 ? true : false, " inner product: ", ortho)
        # println("is Rj*Rj' positive definite? Eigenvalues: ", eigvals(Rj'*Rj))
        # Answer: yes!!!
        # Then we don't need to do QR! We basically just have to normalize the columns of Rj, because off diagonal elements in the R matrix are proportional to the inner product of columns of Rj
        F = qr(Rj)
        # new Q
        # S = diagm(sign.(real.(diag(F.R))))
        # S = I
        # we have to choose S s.t. the symmetry properties are satisfied
        S = sign(real(F.R[1, 1]))
        # Ss[j, :, :] = S
        R = S*F.R
        newQj = Matrix(F.Q)*S
        # println(typeof(newQj))
        # # println(size(Rj))
        # Rjnorms = zeros(Complex{BigFloat}, size(Rj, 2))
        # for (i, col) in enumerate(eachcol(Rj))
        #     colnorm = norm(col)
        #     Rjnorms[i] = colnorm
        # end
        # R = diagm(Rjnorms)
        # foreach(normalize!, eachcol(Rj))
        # newQj = Rj
        # reorthogonalize
        newQj = ChannelMixingDiscretization.gramschmidt(newQj, Q[:, n+1:end])
        # println("new Q block vector orthogonalization: ", newQj'*newQj)
        Q = hcat(Q, newQj)
        
        if (j != m) && norm(R, 1) < 1e-20
            @warn "WARNING! Bad Krylov space!"
        end
        
        # store the results
        Al[j, :, :] = Aj
        B[j + 1, :, :] = R
        Bu[j, :, :] = R'
    end
    Bl = B[2:end-1, :, :]
    # return a dictionary of arrays; keys denote the (off)-diagonal index
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end

function tridiagonalizeSPSU2(A, q, m)
    n = size(q, 2)
    # set the number of iterations as 2J
    m = m === nothing ? Int(size(A)[1]/n) : m
    # Ss = zeros(Float64, m, n, n)
    Al = zeros(Complex{BigFloat}, m, n, n)
    B = zeros(Complex{BigFloat}, m+1, n, n)
    Bu = zeros(Complex{BigFloat}, m, n, n)
    Q = hcat(zeros(BigFloat, size(q)), q)
    for j in 1:m
        Qj = Q[:, end-(n-1):end]
        Uj = A*Qj - Q[:, (end - 2n + 1):(end - n)]*(B[j, :, :]')
        Aj = Qj'*Uj
        Rj = Uj - Qj*Aj
        # ortho = dot(Rj[:, 1], Rj[:, 2])
        # println("are the columns of Rj orthogonal? ", norm(ortho) < 1e-20 ? true : false, " inner product: ", ortho)
        # println("is Rj*Rj' positive definite? Eigenvalues: ", eigvals(Rj'*Rj))
        # Answer: yes!!!
        # Then we don't need to do QR! We basically just have to normalize the columns of Rj, because off diagonal elements in the R matrix are proportional to the inner product of columns of Rj
        F = qr(Rj)
        # new Q
        # S = diagm(sign.(real.(diag(F.R))))
        # S = I
        # we have to choose S s.t. the symmetry properties are satisfied
        S = diagm([sign(real(F.R[1, 1])) < 0 ? -1. : 1., sign(real(F.R[2, 2])) < 0 ? 1. : -1.])
        # Ss[j, :, :] = S
        R = S*F.R
        newQj = Matrix(F.Q)*S
        # println(typeof(newQj))
        # # println(size(Rj))
        # Rjnorms = zeros(Complex{BigFloat}, size(Rj, 2))
        # for (i, col) in enumerate(eachcol(Rj))
        #     colnorm = norm(col)
        #     Rjnorms[i] = colnorm
        # end
        # R = diagm(Rjnorms)
        # foreach(normalize!, eachcol(Rj))
        # newQj = Rj
        # reorthogonalize
        newQj = ChannelMixingDiscretization.gramschmidt(newQj, Q[:, n+1:end])
        # println("new Q block vector orthogonalization: ", newQj'*newQj)
        Q = hcat(Q, newQj)
        
        if (j != m) && norm(R, 1) < 1e-20
            println("WARNING! Bad Krylov space!")
        end
        
        # store the results
        Al[j, :, :] = Aj
        B[j + 1, :, :] = R
        Bu[j, :, :] = R'
    end
    Bl = B[2:end-1, :, :]
    # return a dictionary of arrays; keys denote the (off)-diagonal index
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end