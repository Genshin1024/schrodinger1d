module Schrodinger1d

using LinearAlgebra:diagm,diag

export solver

"""
    solver(potential,n=1, mass=1, grids=300, xlims=[-4,4])

Solve 1D schrodinger equation by finite difference method.
"""
function solver(potential, n=1, mass=1, grids=300,xlims=[-4,4])
    hbar = 1
    H = zeros(grids,grids)
    
    x = range(xlims[1],xlims[2],length=grids)
    U = diagm(map(potential,x)) 
    lap = ( diagm(1 => ones(grids-1), -1 => ones(grids-1)) + -2*diagm(ones(grids)) ) / ((x[2] -x[1])^2)
    H = -hbar^2/(2*mass) * lap + U

    v, d = eigs_qr(H)
    en = grids-n+1
    energy = d[en]
    psi = v[:, en]

    return energy, psi
    
end

"""
    eigs_qr(H)

Solve the eigenvaues problem by QR iteration method.
"""
function eigs_qr(H)
    iters = 300
    tol = 1e-11
    n,m = size(H)
    v = zeros(n,m)
    d = zeros(n)
    s = diagm(ones(n))

    for _ in 1:iters
        q, r = my_qr(H)
        H0 = r*q
        s1 = s*q

        if abs(maximum(diag(H0)) - maximum(diag(H))) < tol
            break
        end

        H = deepcopy(H0)
        s = deepcopy(s1)
    end
    d = diag(H)
    v = s
    
    return v, d
end

"""
    my_qr(H)

Compute the QR factorization of the matrix `H`.
"""
function my_qr(H)
    tol = 1e-7
    r = H
    grids= size(H, 1)
    q = diagm(ones(grids))

    for j in 1:grids 
        for i in j+1:grids
            if abs(r[i,j] -0) < tol
                continue
            end

            a = r[j,j]
            b = r[i,j]

            base = sqrt(a^2+b^2)
            c = a/base
            s = b/base

            g = diagm(ones(grids))
            g[j,j] = c
            g[i,j] = -s
            g[j,i] = s
            g[i,i] = c
            
            r = deepcopy(g*r)
            q = deepcopy(g*q)
        end 
    end
    q = q'

    return q,r
end

end # module
