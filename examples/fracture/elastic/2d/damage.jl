function build_phase_field_sys(H, p::Params)

    Nx, Ny = p.Nx, p.Ny
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
    # for i=1; forward FD
    for j in 2:Ny-1
        k = lin(1,j,Nx)
        M[k,k] = -1/Δx
        M[k,lin(2,j,Nx)] = 1/Δx
    end
    # for i=Nx; backward FD
    for j in 2:Ny-1
        k = lin(Nx,j,Nx)
        M[k,k] = 1/Δx
        M[k,lin(Nx-1,j,Nx)] = -1/Δx
    end
    
    # corner i=1,j=1
    k = lin(1,1,Nx)
    #M[k,k] += -1/Δx
    #M[k,lin(2,1,Nx)] += 1/Δx
    M[k,k] += -1/Δy
    M[k,lin(1,2,Nx)] += 1/Δy
    
    # corner i=1,j=Ny
    k = lin(1,Ny,Nx)
    M[k,k] += -1/Δx
    M[k,lin(2,Ny,Nx)] += 1/Δx
    # M[k,k] += 1/Δy
    # M[k,lin(1,Ny-1,Nx)] += -1/Δy

    # corner i=Nx,j=1
    k = lin(Nx,1,Nx)
    M[k,k] += 1/Δx
    M[k,lin(Nx-1,1,Nx)] += -1/Δx
    # M[k,k] += -1/Δy
    # M[k,lin(Nx,2,Nx)] += 1/Δy

    # corner i=Nx,j=Ny
    k = lin(Nx,Ny,Nx)
    M[k,k] += 1/Δx
    M[k,lin(Nx-1,Ny,Nx)] += -1/Δx
    # M[k,k] += 1/Δy
    # M[k,lin(Nx,Ny-1,Nx)] += -1/Δy

    return M, b
end

#########################
# from chatgpt
function build_phase_field_sys_auto(H, p::Params)

    Nx, Ny = p.Nx, p.Ny
    Δx, Δy = p.Δx, p.Δy
    l = p.l
    C3 = p.C3

    Ntot = Nx * Ny

    M = zeros(Ntot, Ntot)
    b = zeros(Ntot)

    # helper
    α(i,j) = 1/l^2 + 2*C3*H[i,j]/l
    rhs(i,j) = -2*C3*H[i,j]/l

    for j in 1:Ny
        for i in 1:Nx

            k = lin(i,j,Nx)

            ax = 1/Δx^2
            ay = 1/Δy^2

            diag = 0.0

            # --- X direction ---
            if i == 1
                # left boundary (Neumann)
                M[k, lin(i+1,j,Nx)] += ax
                diag -= ax
            elseif i == Nx
                # right boundary (Neumann)
                M[k, lin(i-1,j,Nx)] += ax
                diag -= ax
            else
                # interior
                M[k, lin(i-1,j,Nx)] += ax
                M[k, lin(i+1,j,Nx)] += ax
                diag -= 2ax
            end

            # --- Y direction ---
            if j == 1
                # bottom boundary (Neumann)
                M[k, lin(i,j+1,Nx)] += ay
                diag -= ay
            elseif j == Ny
                # top boundary (Neumann)
                M[k, lin(i,j-1,Nx)] += ay
                diag -= ay
            else
                # interior
                M[k, lin(i,j-1,Nx)] += ay
                M[k, lin(i,j+1,Nx)] += ay
                diag -= 2ay
            end

            # --- reaction term ---
            diag -= α(i,j)

            M[k,k] += diag
            b[k] = rhs(i,j)

        end
    end

    return M, b
end
#########################

# copy to the ghost zone the second to last point, s.t. the centred FC is zero, both x,y, that is f_0 = f_2 -> df_1 = (f_2 - f_0)/ h
function build_phase_field_sys_ghost(H, p::Params)

    Nx, Ny = p.Nx, p.Ny
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
    
    # bottom boundary
    j=1
    for i in 2:Nx-1
        k = lin(i,j,Nx)
        M[k,lin(i-1,j,Nx)] = 1/Δx^2
        #M[k,lin(i,j-1,Nx)] = 1/Δy^2
        M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δy^2
        M[k,lin(i,j+1,Nx)] = 1/Δy^2
        M[k,lin(i+1,j,Nx)] = 1/Δx^2
        b[k] = -2*C3*H[i,j]/l
    end

    # top boundary
    j=Ny
    for i in 2:Nx-1
        k = lin(i,j,Nx)
        M[k,lin(i-1,j,Nx)] = 1/Δx^2
        M[k,lin(i,j-1,Nx)] = 1/Δy^2
        M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δy^2
        #M[k,lin(i,j+1,Nx)] = 1/Δy^2
        M[k,lin(i+1,j,Nx)] = 1/Δx^2
        b[k] = -2*C3*H[i,j]/l
    end

    # left boundary
    i=1
    for j in 2:Ny-1
        k = lin(i,j,Nx)
        #M[k,lin(i-1,j,Nx)] = 1/Δx^2
        M[k,lin(i,j-1,Nx)] = 1/Δy^2
        M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δx^2
        M[k,lin(i,j+1,Nx)] = 1/Δy^2
        M[k,lin(i+1,j,Nx)] = 1/Δx^2
        b[k] = -2*C3*H[i,j]/l
    end

    # right boundary
    i=Nx
    for j in 2:Ny-1
        k = lin(i,j,Nx)
        M[k,lin(i-1,j,Nx)] = 1/Δx^2
        M[k,lin(i,j-1,Nx)] = 1/Δy^2
        M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δx^2
        M[k,lin(i,j+1,Nx)] = 1/Δy^2
        #M[k,lin(i+1,j,Nx)] = 1/Δx^2
        b[k] = -2*C3*H[i,j]/l
    end

    #corner
    i = 1
    j = 1
    k = lin(i,j,Nx)
    #M[k,lin(i-1,j,Nx)] = 1/Δx^2
    #M[k,lin(i,j-1,Nx)] = 1/Δy^2
    M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δx^2 + 1/Δy^2
    M[k,lin(i,j+1,Nx)] = 1/Δy^2
    M[k,lin(i+1,j,Nx)] = 1/Δx^2
    b[k] = -2*C3*H[i,j]/l

    #corner
    i = 1
    j = Ny
    k = lin(i,j,Nx)
    #M[k,lin(i-1,j,Nx)] = 1/Δx^2
    M[k,lin(i,j-1,Nx)] = 1/Δy^2
    M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δx^2 + 1/Δy^2
    #M[k,lin(i,j+1,Nx)] = 1/Δy^2
    M[k,lin(i+1,j,Nx)] = 1/Δx^2
    b[k] = -2*C3*H[i,j]/l

    #corner
    i = Nx
    j = 1
    k = lin(i,j,Nx)
    M[k,lin(i-1,j,Nx)] = 1/Δx^2
    #M[k,lin(i,j-1,Nx)] = 1/Δy^2
    M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δy^2 + 1/Δx^2
    M[k,lin(i,j+1,Nx)] = 1/Δy^2
    #M[k,lin(i+1,j,Nx)] = 1/Δx^2
    b[k] = -2*C3*H[i,j]/l

    #corner
    i = Nx
    j = Ny
    k = lin(i,j,Nx)
    M[k,lin(i-1,j,Nx)] = 1/Δx^2
    M[k,lin(i,j-1,Nx)] = 1/Δy^2
    M[k,k] = -2/Δx^2 - 2/Δy^2 - 1/l^2 - 2*C3*H[i,j]/l + 1/Δy^2 + 1/Δx^2
    #M[k,lin(i,j+1,Nx)] = 1/Δy^2
    #M[k,lin(i+1,j,Nx)] = 1/Δx^2
    b[k] = -2*C3*H[i,j]/l

    return M, b
end

function history_from_d(d, p::Params,f1_left,f1_right)

    Nx, Ny = p.Nx, p.Ny
    H = zeros(Nx,Ny)

    for j in 1:Ny
        for i in 1:Nx
            H[i,j] = (d[i,j] -p.l^2*(Dxx_ghost(d,i,j,p,f1_left,f1_right) + Dyy_ghost(d,i,j,p)))/(2*p.l*p.C3*(1-d[i,j]+p.k))
        end
    end

    return H
end

# function history_from_d(d, p::Params)

#     Nx, Ny = p.Nx, p.Ny
#     H = zeros(Nx,Ny)

#     for j in 2:Ny-1
#         for i in 2:Nx-1
#             H[i,j] = (d[i,j] -p.l^2*(Dxx(d,i,j,p) + Dyy(d,i,j,p)))/(2*p.l*p.C3*(1-d[i,j]+p.k))
#         end
#     end

#     j = 1
#     for i in 1:Nx
#         H[i,j] = (d[i,j])/(2*p.l*p.C3*(1-d[i,j]+p.k))
#     end

#     j = Ny
#     for i in 1:Nx
#         H[i,j] = (d[i,j])/(2*p.l*p.C3*(1-d[i,j]+p.k))
#     end

#     i = 1
#     for j in 2:Ny-1
#         H[i,j] = (d[i,j])/(2*p.l*p.C3*(1-d[i,j]+p.k))
#     end

#     i = Nx
#     for j in 2:Ny-1
#         H[i,j] = (d[i,j])/(2*p.l*p.C3*(1-d[i,j]+p.k))
#     end

#     return H
# end

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

# knows of ghost zones
function history_ghost(H_old, a1, a2, p::Params,f1_left,f1_right)
    λ, μ = p.λ, p.μ
    H_new = zeros(p.Nx,p.Ny)
    for j in 1:p.Ny
        for i in 1:p.Nx
            H_new[i,j] = (λ + 2*μ/2)*(Dx_ghost(a1,i,j,p,f1_left,f1_right) + Dy_ghost(a2,i,j,p))^2
        end
    end
    return max.(H_old, H_new)
end