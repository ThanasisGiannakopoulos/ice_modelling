function build_phase_field_sys(H, p::Params)

    Nx, Ny = p.Nx, p.Nx
    Δx, Δy = p.Δx, p.Δy
    l = p.l
    C3 = p.C3
    
    Ntot = Nx*Ny # d vectorized

    M = zeros(Ntot,Ntot)
    b = zeros(Ntot)

    for j in 2:Ny-1
        for i in 2:Nx-1
            k = lin(i,j,Nx)
            M[k,lin(i-1,j,Nx)] = 1/Δx^2
            M[k,lin(i,j-1,Nx)] = 1/Δy^2
            M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l
            M[k,lin(i,j+1,Nx)] = 1/Δy^2
            M[k,lin(i+1,j,Nx)] = 1/Δx^2
            b[k] = -2*C3*H[i,j]/l
        end
    end
    # BC is d*n = 0 -> d_x=0 on left/right and d_y=0 on top/bottom
    # Neumann with back/forward 1st order FD
    # for j=1; forward FD
    for i in 2:Nx-1
        k = lin(i,1,Nx)
        M[k,k] = -1/Δy
        M[k,lin(i,2,Nx)] = 1/Δy
    end
    # for j=Ny; backward FD
    for i in 2:Nx-1
        k = lin(i,Ny,Nx)
        M[k,k] = 1/Δy
        M[k,lin(i,Ny-1,Nx)] = -1/Δy
    end
    # for i=1; forward FD; includes corners i=1,j=1 and i=1,j=Ny
    for j in 1:Ny
        k = lin(1,j,Nx)
        M[k,k] = -1/Δx
        M[k,lin(2,j,Nx)] = 1/Δx
    end
    # for i=Nx; backward FD; includes corners i=Nx,j=1 and i=Nx,j=Ny
    for j in 1:Ny
        k = lin(Nx,j,Nx)
        M[k,k] = 1/Δx
        M[k,lin(Nx-1,j,Nx)] = -1/Δx
    end
    
    return M, b
end

function history_from_d(d, p::Params)

    Nx, Ny = p.Nx, p.Ny
    H = zeros(Nx,Ny)

    for j in 1:Ny
        for i in 1:Nx
            H[i,j] = (d[i,j] -p.l^2*(Dxx(d,i,j,p) + Dyy(d,i,j,p)))/(2*p.l*p.C3*(1-d[i,j]+p.k))
        end
    end

    return H
end

# just (λ + 2*μ/2)*tr(e)^2 for now, no splits etc, just to test 
# probably ok for just pull (tensile) in x
function history(H_old, a1, a2, p::Params)
    λ, μ = p.λ, p.μ
    H_new = zeros(p.Nx,p.Ny)
    for j in 1:p.Ny
        for i in 1:p.Nx
            H_new[i,j] = (λ + 2*μ/2)*(Dx(a1,i,j,p) + Dy(a2,i,j,p))^2
        end
    end
    return max.(H_old, H_new)
end