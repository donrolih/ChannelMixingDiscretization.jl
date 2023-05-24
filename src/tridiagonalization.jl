function gramschmidt(u, Q; maxiter=3, α=big(0.5))
    r = norm(u'*u)
    r1 = 0.
    u_start = u
    for i in 1:maxiter
        u = u - Q*Q'*u
        r1 = norm(u'*u)
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
    # return a dictionary of arrays; keys denote the (off)diagonal 
    diags = Dict(-1 => Bl, 0 => Al, 1 => Bu[1:end-1, :, :])
    return diags
end
