# Build global matrix M and rhs b
function build_displacement_system(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p::Params,
    f1_left, f1_right,
    f2_left, f2_right)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny # a1, a2 are vectorized
    Ndof = 2*Ntot # for both a1 and a2

    λ, μ = p.λ, p.μ

    M = zeros(Ndof,Ndof)
    b = zeros(Ndof)

    # BC in x
    # left
    if p.left_x_BC=="Dirichlet"
        # for i = 1
        for j in 1:Ny
            k = lin(1,j,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            M[k1,k1]=1; b[k1]=f1_left[j] # grid function of y
            M[k2,k2]=1; b[k2]=f2_left[j]
        end
    elseif p.left_x_BC=="traction"
        # for i = 1
        # for j=2 we apply another constraint to remove rotation
        for j in 2:Ny-1
            k = lin(1,j,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = f1_left
            # tensile/compressive
            # (λ+2μ)∂_x a1 # forward FD
            M[k1,k1]= (λ+2μ)*(-1/Δx)
            M[k1,lin(2,j,Nx)]= (λ+2μ)*(1/Δx) 
            # + λ ∂_y a2 # centred FD
            M[k1,Ntot+lin(1,j-1,Nx)] = λ*(-0.5/Δy) #μ*(-0.5/Δy)
            M[k1,Ntot+lin(1,j+1,Nx)] = λ*(0.5/Δy) # μ*(0.5/Δy)
            # prescribed f1_left
            b[k1] = -f1_left[j] 
            # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
            # μ ∂_y a1; centred FD
            M[k2,lin(1,j-1,Nx)] = μ*(-0.5/Δy)
            M[k2,lin(1,j+1,Nx)] = μ*(0.5/Δy)
            # + μ ∂_x a2; forward FD
            M[k2,Ntot+lin(1,j,Nx)] = μ*(-1/Δx)
            M[k2,Ntot+lin(2,j,Nx)] = μ*(1/Δx)
            # prescribed f2_left
            b[k2] = -f2_left[j]
        end
        
        # -------------------------------------------------------------------------
        # on the 2 corners use forward and backward FD for ∂_y
        # -------------------------------------------------------------------------
        # i=1,j=1
        k = lin(1,1,Nx)
        k1 = k # for a1 positions in vector [a1, a2]
        k2 = k + Ntot # for a2...
        # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = f1_left
        # tensile/compressive
        # (λ+2μ)∂_x a1 # forward FD
        M[k1,k1]= (λ+2μ)*(-1/Δx)
        M[k1,lin(2,1,Nx)]= (λ+2μ)*(1/Δx) 
        # + λ ∂_y a2 # forward FD
        M[k1,Ntot+lin(1,1,Nx)] = λ*(-1/Δy) #μ*(-1/Δy)
        M[k1,Ntot+lin(1,2,Nx)] = λ*(1/Δy) #μ*(1/Δy)
        # prescribed f1_left
        b[k1] = -f1_left[1] 
        # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
        # μ ∂_y a1; forward FD
        M[k2,lin(1,1,Nx)] = μ*(-1/Δy)
        M[k2,lin(1,2,Nx)] = μ*(1/Δy)
        # + μ ∂_x a2; forward FD
        M[k2,Ntot+lin(1,1,Nx)] = μ*(-1/Δx)
        M[k2,Ntot+lin(2,1,Nx)] = μ*(1/Δx)
        # prescribed f2_left
        b[k2] = -f2_left[1]
        
        # -------------------------------------------------------------------------
        # i=1,j=Ny
        # -------------------------------------------------------------------------
        
        k = lin(1,Ny,Nx)
        k1 = k # for a1 positions in vector [a1, a2]
        k2 = k + Ntot # for a2...
        # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = f1_left
        # tensile/compressive
        # (λ+2μ)∂_x a1 # forward FD
        M[k1,k1]= (λ+2μ)*(-1/Δx)
        M[k1,lin(2,Ny,Nx)]= (λ+2μ)*(1/Δx) 
        # + λ ∂_y a2 # backward FD
        M[k1,Ntot+lin(1,Ny-1,Nx)] = λ*(-1/Δy) #μ*(-1/Δy)
        M[k1,Ntot+lin(1,Ny,Nx)] = λ*(1/Δy) #μ*(1/Δy)
        # prescribed f1_left
        b[k1] = -f1_left[Ny] 
        # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
        # μ ∂_y a1; backward FD
        M[k2,lin(1,Ny-1,Nx)] = μ*(-1/Δy)
        M[k2,lin(1,Ny,Nx)] = μ*(1/Δy)
        # + μ ∂_x a2; forward FD
        M[k2,Ntot+lin(1,Ny,Nx)] = μ*(-1/Δx)
        M[k2,Ntot+lin(2,Ny,Nx)] = μ*(1/Δx)
        # prescribed f2_left
        b[k2] = -f2_left[Ny]
    
    else
        println("no such left_x_BC")
    end
    # right BC
    if p.right_x_BC=="Dirichlet"
        # for i = Nx
        for j in 1:Ny
            k = lin(Nx,j,Nx)
            k1 = k # for a1
            k2 = k + Ntot # for a2
            M[k1,k1]=1; b[k1]=f1_right[j] # grid function of y
            M[k2,k2]=1; b[k2]=f2_right[j]
        end
    elseif p.right_x_BC=="traction"
        # for i = Nx
        for j in 2:Ny-1
            k = lin(Nx,j,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = t11
            # tensile/compressive
            # (λ+2μ)∂_x a1 # backward FD
            M[k1,k1]= (λ+2μ)*(1/Δx);
            M[k1,lin(Nx-1,j,Nx)]= (λ+2μ)*(-1/Δx); 
            # + λ ∂_y a2 # centred FD
            M[k1,Ntot+lin(Nx,j-1,Nx)] = λ*(-0.5/Δy) # μ*(-0.5/Δy)
            M[k1,Ntot+lin(Nx,j+1,Nx)] = λ*(0.5/Δy) # μ*(0.5/Δy)
            # prescribed f1_right
            b[k1] = f1_right[j] 
            # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
            # μ ∂_y a1; centred FD
            M[k2,lin(Nx,j-1,Nx)] = μ*(-0.5/Δy)
            M[k2,lin(Nx,j+1,Nx)] = μ*(0.5/Δy)
            # + μ ∂_x a2; backward FD
            M[k2,Ntot+lin(Nx-1,j,Nx)] = μ*(-1/Δx)
            M[k2,Ntot+lin(Nx,j,Nx)] = μ*(1/Δx)
            # prescribed f2_right
            b[k2] = f2_right[j]
        end
        # for the corners
        # i=Nx,j=1
        k = lin(Nx,1,Nx)
        k1 = k # for a1 positions in vector [a1, a2]
        k2 = k + Ntot # for a2...
        # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = t11
        # tensile/compressive
        # (λ+2μ)∂_x a1 # backward FD
        M[k1,k1]= (λ+2μ)*(1/Δx);
        M[k1,lin(Nx-1,1,Nx)]= (λ+2μ)*(-1/Δx); 
        # + λ ∂_y a2 # forward FD
        M[k1,Ntot+lin(Nx,1,Nx)] = λ*(-1/Δy) #μ*(-1/Δy)
        M[k1,Ntot+lin(Nx,2,Nx)] = λ*(1/Δy) #μ*(1/Δy)
        # prescribed f1_right
        b[k1] = f1_right[1] 
        # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
        # μ ∂_y a1; forward FD
        M[k2,lin(Nx,1,Nx)] = μ*(-1/Δy)
        M[k2,lin(Nx,2,Nx)] = μ*(1/Δy)
        # + μ ∂_x a2; backward FD
        M[k2,Ntot+lin(Nx-1,1,Nx)] = μ*(-1/Δx)
        M[k2,Ntot+lin(Nx,1,Nx)] = μ*(1/Δx)
        # prescribed f2_right
        b[k2] = f2_right[1]
        # i=Nx,j=Ny
        k = lin(Nx,Ny,Nx)
        k1 = k # for a1 positions in vector [a1, a2]
        k2 = k + Ntot # for a2...
        # for sigma_11 = (λ+2μ)∂_x a1 + λ ∂_y a2 = t11
        # tensile/compressive
        # (λ+2μ)∂_x a1 # backward FD
        M[k1,k1]= (λ+2μ)*(1/Δx);
        M[k1,lin(Nx-1,Ny,Nx)]= (λ+2μ)*(-1/Δx); 
        # + λ ∂_y a2 # backward FD
        M[k1,Ntot+lin(Nx,Ny-1,Nx)] = λ*(-1/Δy) # μ*(-1/Δy)
        M[k1,Ntot+lin(Nx,Ny,Nx)] = λ*(1/Δy) # μ*(1/Δy)
        # prescribed f1_right
        b[k1] = f1_right[Ny] 
        # for sigma_12 = μ ∂_y a1 + μ ∂_x a2 = f2_left
        # μ ∂_y a1; backward FD
        M[k2,lin(Nx,Ny-1,Nx)] = μ*(-1/Δy)
        M[k2,lin(Nx,Ny,Nx)] = μ*(1/Δy)
        # + μ ∂_x a2; backward FD
        M[k2,Ntot+lin(Nx-1,Ny,Nx)] = μ*(-1/Δx)
        M[k2,Ntot+lin(Nx,Ny,Nx)] = μ*(1/Δx)
        # prescribed f2_right
        b[k2] = f2_right[Ny]
    else
        println("no such right_x_BC")
    end

    
    # BC in y
    if p.y_BC=="Neumann" # using 1st order stencils
        # for j=1
        for i in 2:Nx-1 # the edges of x are given by the Dirichlet BCs in x
            k = lin(i,1,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            M[k1,k1] = -1/Δy
            M[k1,lin(i,2,Nx)] = 1/Δy
            M[k2,k2] = -1/Δy
            M[k2,lin(i,2,Nx)+Ntot] = 1/Δy
        end
        # for j=Ny
        for i in 2:Nx-1 # the edges of x are given by the Dirichlet BCs in x
            k = lin(i,Ny,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            M[k1,k1] = 1/Δy
            M[k1,lin(i,Ny-1,Nx)] = -1/Δy
            M[k2,k2] = 1/Δy
            M[k2,lin(i,Ny-1,Nx)+Ntot] = -1/Δy
        end
    elseif p.y_BC == "traction_free"
        # for j=1
        for i in 2:Nx-1 # the edges of x are given by the Dirichlet BCs in x
            k = lin(i,1,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            # for sigma_12 = 0 (μ*∂_y*a1 + μ*∂_x*a2 = 0)
            M[k1,k1] = μ*(-1/Δy)
            M[k1,lin(i,2,Nx)] = μ*(1/Δy)
            M[k1,lin(i+1,1,Nx)+Ntot] = μ*(1/(2*Δx))
            M[k1,lin(i-1,1,Nx)+Ntot] = μ*(-1/(2*Δx))
            # for sigma_22 = 0 (λ*∂_x*a1 + (λ+2μ)*∂_y*a2 = 0)
            M[k2,lin(i+1,1,Nx)] = λ*(1/(2*Δx))
            M[k2,lin(i-1,1,Nx)] = λ*(-1/(2*Δx))
            M[k2,k2] = (λ+2*μ)*(-1/Δy)
            M[k2,lin(i,2,Nx)+Ntot] = (λ+2*μ)*(1/Δy)
        end
        # for j=Ny
        for i in 2:Nx-1 # the edges of x are given by the Dirichlet BCs in x
            k = lin(i,Ny,Nx)
            k1 = k # for a1 positions in vector [a1, a2]
            k2 = k + Ntot # for a2...
            # for sigma_12 = 0 (μ*∂_y*a1 + μ*∂_x*a2 = 0)
            M[k1,k1] = μ*(1/Δy)
            M[k1,lin(i,Ny-1,Nx)] = μ*(-1/Δy)
            M[k1,lin(i+1,Ny,Nx)+Ntot] = μ*(1/(2*Δx))
            M[k1,lin(i-1,Ny,Nx)+Ntot] = μ*(-1/(2*Δx))
            # for sigma_22 = 0 (λ*∂_x*a1 + (λ+2μ)*∂_y*a2 = 0)
            M[k2,lin(i+1,Ny,Nx)] = λ*(1/(2*Δx))
            M[k2,lin(i-1,Ny,Nx)] = λ*(-1/(2*Δx))
            M[k2,k2] = (λ+2*μ)*(1/Δy)
            M[k2,lin(i,Ny-1,Nx)+Ntot] = (λ+2*μ)*(-1/Δy)

        end
    else
        println("no such y_BC")
    end
    
    for j in 2:Ny-1
        for i in 2:Nx-1

            # midpoints in x
            # sigma_11
            P111_i_mhalf_j = 0.5*(P111[i-1,j] + P111[i,j])
            P111_i_phalf_j = 0.5*(P111[i,j]   + P111[i+1,j])
            Q111_i_mhalf_j = 0.5*(Q111[i-1,j] + Q111[i,j])
            Q111_i_phalf_j = 0.5*(Q111[i,j]   + Q111[i+1,j])
            P112_i_mhalf_j = 0.5*(P112[i-1,j] + P112[i,j])
            P112_i_phalf_j = 0.5*(P112[i,j]   + P112[i+1,j])
            Q112_i_mhalf_j = 0.5*(Q112[i-1,j] + Q112[i,j])
            Q112_i_phalf_j = 0.5*(Q112[i,j]   + Q112[i+1,j])
            # sigma_12
            P121_i_mhalf_j = 0.5*(P121[i-1,j] + P121[i,j])
            P121_i_phalf_j = 0.5*(P121[i,j]   + P121[i+1,j])
            Q121_i_mhalf_j = 0.5*(Q121[i-1,j] + Q121[i,j])
            Q121_i_phalf_j = 0.5*(Q121[i,j]   + Q121[i+1,j])
            P122_i_mhalf_j = 0.5*(P122[i-1,j] + P122[i,j])
            P122_i_phalf_j = 0.5*(P122[i,j]   + P122[i+1,j])
            Q122_i_mhalf_j = 0.5*(Q122[i-1,j] + Q122[i,j])
            Q122_i_phalf_j = 0.5*(Q122[i,j]   + Q122[i+1,j])

            # midpoints in y
            # sigma_12
            P121_i_j_mhalf = 0.5*(P121[i,j-1] + P121[i,j])
            P121_i_j_phalf = 0.5*(P121[i,j]   + P121[i,j+1])
            Q121_i_j_mhalf = 0.5*(Q121[i,j-1] + Q121[i,j])
            Q121_i_j_phalf = 0.5*(Q121[i,j]   + Q121[i,j+1])
            P122_i_j_mhalf = 0.5*(P122[i,j-1] + P122[i,j])
            P122_i_j_phalf = 0.5*(P122[i,j]   + P122[i,j+1])
            Q122_i_j_mhalf = 0.5*(Q122[i,j-1] + Q122[i,j])
            Q122_i_j_phalf = 0.5*(Q122[i,j]   + Q122[i,j+1])
            # sigma_22
            P221_i_j_mhalf = 0.5*(P221[i,j-1] + P221[i,j])
            P221_i_j_phalf = 0.5*(P221[i,j]   + P221[i,j+1])
            Q221_i_j_mhalf = 0.5*(Q221[i,j-1] + Q221[i,j])
            Q221_i_j_phalf = 0.5*(Q221[i,j]   + Q221[i,j+1])
            P222_i_j_mhalf = 0.5*(P222[i,j-1] + P222[i,j])
            P222_i_j_phalf = 0.5*(P222[i,j]   + P222[i,j+1])
            Q222_i_j_mhalf = 0.5*(Q222[i,j-1] + Q222[i,j])
            Q222_i_j_phalf = 0.5*(Q222[i,j]   + Q222[i,j+1])

            k = lin(i,j,Nx)
            k1 = k # for a1
            k2 = k + Ntot # for a2

            # 1st eq; rhs = 0
            # for a1 part
            # a1_i-1_j-1
            M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
            # a1_i-1_j
            M[k1,lin(i-1,j,Nx)] = Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
            # a1_i-1_j+1
            M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
            # a1_i_j-1
            M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
            # a1_i,j
            M[k1,k1] = -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
            # a1_i_j+1
            M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
            # a1_i+1_j-1
            M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
            # a1_i+1,j
            M[k1,lin(i+1,j,Nx)] = Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
            # a1_i+1_j+1
            M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
            
            # for a2 part
            # a2_i-1_j-1
            M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
            # a2_i-1_j
            M[k1,lin(i-1,j,Nx)+Ntot] = Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
            # a2_i-1_j+1
            M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
            # a2_i_j-1
            M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
            # a2_i_j
            M[k1,k2] = -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
            # a2_i_j+1
            M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
            # a2_i+1_j-1
            M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
            # a2_i+1_j
            M[k1,lin(i+1,j,Nx)+Ntot] = Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
            # a2_i+1_j+1
            M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)

            # 2nd eq; rhs = 0
            # for a1 part
            # a1_i-1_j-1
            M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
            # a1_i-1_j
            M[k2,lin(i-1,j,Nx)] = Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
            # a1_i-1_j+1
            M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
            # a1_i_j-1
            M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
            # a1_i,j
            M[k2,k1] = -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
            # a1_i_j+1
            M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
            # a1_i+1_j-1
            M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
            # a1_i+1,j
            M[k2,lin(i+1,j,Nx)] = Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
            # a1_i+1_j+1
            M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
            
            # for a2 part
            # a2_i-1_j-1
            M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
            # a2_i-1_j
            M[k2,lin(i-1,j,Nx)+Ntot] = Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
            # a2_i-1_j+1
            M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
            # a2_i_j-1
            M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
            # a2_i_j
            M[k2,k2] = -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
            # a2_i_j+1
            M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
            # a2_i+1_j-1
            M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
            # a2_i+1_j
            M[k2,lin(i+1,j,Nx)+Ntot] = Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
            # a2_i+1_j+1
            M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)

        end
    end

    return M,b
end

function displacement_BC_left_right(t::Float64, p::Params)
    f1_left  = -t*1e-2*ones(p.Ny) # 0.0*ones(p.Ny)#-t*1.0*1e0*ones(p.Ny) # 0.0*ones(p.Ny)
    f2_left  = 0.0*ones(p.Ny)
    
    f1_right =  0.0*t*0.5*1e-2*ones(p.Ny)
    f2_right =  0*1e-1*ones(p.Ny)
    return f1_left, f2_left, f1_right, f2_right
end