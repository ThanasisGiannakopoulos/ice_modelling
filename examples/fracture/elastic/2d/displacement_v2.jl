# Build global matrix M and rhs b
function build_displacement_system(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p::Params,
    f1_left, f1_right,
    f2_left, f2_right)

    Nx,Ny = p.Nx ,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny # a1, a2 are vectorized  # remove ghost zones from left-right and top-bottom
    Ndof = 2*Ntot # for both a1 and a2

    λ, μ = p.λ, p.μ

    M = zeros(Ndof,Ndof)
    b = zeros(Ndof)
    
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

    # i=1, j=2:Ny-1 i.e. left boundary, no corners; ghost points on the left for Dirichlet BC
    for j in 2:Ny-1
        i = 1
        # midpoints in x
        # sigma_11
        P111_i_mhalf_j = P111[i,j]
        P111_i_phalf_j = 0.5*(P111[i,j] + P111[i+1,j])
        Q111_i_mhalf_j = Q111[i,j]
        Q111_i_phalf_j = 0.5*(Q111[i,j] + Q111[i+1,j])
        P112_i_mhalf_j = P112[i,j]
        P112_i_phalf_j = 0.5*(P112[i,j] + P112[i+1,j])
        Q112_i_mhalf_j = Q112[i,j]
        Q112_i_phalf_j = 0.5*(Q112[i,j] + Q112[i+1,j])
        # sigma_12
        P121_i_mhalf_j = P121[i,j]
        P121_i_phalf_j = 0.5*(P121[i,j] + P121[i+1,j])
        Q121_i_mhalf_j = Q121[i,j]
        Q121_i_phalf_j = 0.5*(Q121[i,j] + Q121[i+1,j])
        P122_i_mhalf_j = P122[i,j]
        P122_i_phalf_j = 0.5*(P122[i,j] + P122[i+1,j])
        Q122_i_mhalf_j = Q122[i,j]
        Q122_i_phalf_j = 0.5*(Q122[i,j] + Q122[i+1,j])

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
        #M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        eq1_a1_im1_jm1 = f1_left[j-1]*(0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy))
        # a1_i-1_j
        #M[k1,lin(i-1,j,Nx)] = Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
        eq1_a1_im1_j = f1_left[j]*(Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy))
        # a1_i-1_j+1
        #M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
        eq1_a1_im1_jp1 = f1_left[j+1]*(-0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy))
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
        #M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        eq1_a2_im1_jm1 = f2_left[j-1]*(0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy))
        # a2_i-1_j
        #M[k1,lin(i-1,j,Nx)+Ntot] = Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
        eq1_a2_im1_j = f2_left[j]*(Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy))
        # a2_i-1_j+1
        #M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
        eq1_a2_im1_jp1 = f2_left[j+1]*(-0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy))
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

        # move replacements to the lhs; pick up a minus
        b[k1] = -eq1_a1_im1_jm1 -eq1_a1_im1_j -eq1_a1_im1_jp1 -eq1_a2_im1_jm1 -eq1_a2_im1_j -eq1_a2_im1_jp1

        # 2nd eq; rhs = 0
        # for a1 part
        # a1_i-1_j-1
        #M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        eq2_a1_im1_jm1 = f1_left[j-1]*(0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy))
        # a1_i-1_j
        #M[k2,lin(i-1,j,Nx)] = Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
        eq2_a1_im1_j = f1_left[j]*(Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy))
        # a1_i-1_j+1
        #M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
        eq2_a1_im1_jp1 = f1_left[j+1]*(-0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy))
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
        #M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        eq2_a2_im1_jm1 = f2_left[j-1]*(0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy))
        # a2_i-1_j
        #M[k2,lin(i-1,j,Nx)+Ntot] = Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
        eq2_a2_im1_j = f2_left[j]*(Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy))
        # a2_i-1_j+1
        #M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
        eq2_a2_im1_jp1 = f2_left[j+1]*(-0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy))
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
        
        # move replacements to the lhs; pick up a minus
        b[k2] = -eq2_a1_im1_jm1 -eq2_a1_im1_j -eq2_a1_im1_jp1 -eq2_a2_im1_jm1 -eq2_a2_im1_j -eq2_a2_im1_jp1

    end # end left Dirichlet BC

    # right Dirichlet BC, no corners
    for j in 2:Ny-1
        i = Nx
        # midpoints in x
        # sigma_11
        P111_i_mhalf_j = 0.5*(P111[i-1,j] + P111[i,j])
        P111_i_phalf_j = P111[i,j]
        Q111_i_mhalf_j = 0.5*(Q111[i-1,j] + Q111[i,j])
        Q111_i_phalf_j = Q111[i,j]
        P112_i_mhalf_j = 0.5*(P112[i-1,j] + P112[i,j])
        P112_i_phalf_j = P112[i,j]
        Q112_i_mhalf_j = 0.5*(Q112[i-1,j] + Q112[i,j])
        Q112_i_phalf_j = Q112[i,j]
        # sigma_12
        P121_i_mhalf_j = 0.5*(P121[i-1,j] + P121[i,j])
        P121_i_phalf_j = P121[i,j]
        Q121_i_mhalf_j = 0.5*(Q121[i-1,j] + Q121[i,j])
        Q121_i_phalf_j = Q121[i,j]
        P122_i_mhalf_j = 0.5*(P122[i-1,j] + P122[i,j])
        P122_i_phalf_j = P122[i,j]
        Q122_i_mhalf_j = 0.5*(Q122[i-1,j] + Q122[i,j])
        Q122_i_phalf_j = Q122[i,j]

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
        #M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        eq1_a1_ip1_jm1 = f1_right[j-1]*(-0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy))
        # a1_i+1,j
        #M[k1,lin(i+1,j,Nx)] = Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
        eq1_a1_ip1_j = f1_right[j]*(Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy))
        # a1_i+1_j+1
        #M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
        eq1_a1_ip1_jp1 = f1_right[j+1]*(0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy))
        
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
        #M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        eq1_a2_ip1_jm1 = f2_right[j-1]*(-0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy))
        # a2_i+1_j
        #M[k1,lin(i+1,j,Nx)+Ntot] = Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
        eq1_a2_ip1_j = f2_right[j]*(Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy))
        # a2_i+1_j+1
        #M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)
        eq1_a2_ip1_jp1 = f2_right[j+1]*(0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy))

        # move replacements to the lhs; pick up a minus
        b[k1] = -eq1_a1_ip1_jm1 -eq1_a1_ip1_j -eq1_a1_ip1_jp1 -eq1_a2_ip1_jm1 -eq1_a2_ip1_j -eq1_a2_ip1_jp1

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
        #M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        eq2_a1_ip1_jm1 = f1_right[j-1]*(-0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy))
        # a1_i+1,j
        #M[k2,lin(i+1,j,Nx)] = Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
        eq2_a1_ip1_j = f1_right[j]*(Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy))
        # a1_i+1_j+1
        #M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
        eq2_a1_ip1_jp1 = f1_right[j+1]*(0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy))
        
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
        #M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        eq2_a2_ip1_jm1 = f2_right[j-1]*(-0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy))
        # a2_i+1_j
        #M[k2,lin(i+1,j,Nx)+Ntot] = Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
        eq2_a2_ip1_j = f2_right[j]*(Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy))
        # a2_i+1_j+1
        #M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
        eq2_a2_ip1_jp1 = f2_right[j+1]*(0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy))
        
        # move replacements to the lhs; pick up a minus
        b[k2] = -eq2_a1_ip1_jm1 -eq2_a1_ip1_j -eq2_a1_ip1_jp1 -eq2_a2_ip1_jm1 -eq2_a2_ip1_j -eq2_a2_ip1_jp1

    end # end right Dirichlet

    # Copy the j=1 points at j=0, i.e. move the coefs of j-1 to j 
    # no corners
    j = 1
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
        P121_i_j_mhalf = P121[i,j]
        P121_i_j_phalf = 0.5*(P121[i,j] + P121[i,j+1])
        Q121_i_j_mhalf = Q121[i,j]
        Q121_i_j_phalf = 0.5*(Q121[i,j] + Q121[i,j+1])
        P122_i_j_mhalf = P122[i,j]
        P122_i_j_phalf = 0.5*(P122[i,j] + P122[i,j+1])
        Q122_i_j_mhalf = Q122[i,j]
        Q122_i_j_phalf = 0.5*(Q122[i,j] + Q122[i,j+1])
        # sigma_22
        P221_i_j_mhalf = P221[i,j]
        P221_i_j_phalf = 0.5*(P221[i,j] + P221[i,j+1])
        Q221_i_j_mhalf = Q221[i,j]
        Q221_i_j_phalf = 0.5*(Q221[i,j] + Q221[i,j+1])
        P222_i_j_mhalf = P222[i,j]
        P222_i_j_phalf = 0.5*(P222[i,j] + P222[i,j+1])
        Q222_i_j_mhalf = Q222[i,j]
        Q222_i_j_phalf = 0.5*(Q222[i,j] + Q222[i,j+1])

        k = lin(i,j,Nx)
        k1 = k # for a1
        k2 = k + Ntot # for a2

        # move all the j-1 coeffs to j (copy j to j-1)
        # 1st eq; rhs = 0
        # for a1 part
        # a1_i-1_j-1
        #M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        M[k1,lin(i-1,j,Nx)] += 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i-1_j
        M[k1,lin(i-1,j,Nx)] += Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
        # a1_i-1_j+1
        M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
        # a1_i_j-1
        #M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
        M[k1,lin(i,j,Nx)] += 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
        # a1_i,j
        M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
        # a1_i_j+1
        M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
        # a1_i+1_j-1
        #M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        M[k1,lin(i+1,j,Nx)] += -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i+1,j
        M[k1,lin(i+1,j,Nx)] += Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i+1_j+1
        M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
        
        # for a2 part
        # a2_i-1_j-1
        #M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        M[k1,lin(i-1,j,Nx)+Ntot] += 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i-1_j
        M[k1,lin(i-1,j,Nx)+Ntot] += Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
        # a2_i-1_j+1
        M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
        # a2_i_j-1
        #M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
        M[k1,lin(i,j,Nx)+Ntot] += 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
        # a2_i_j
        M[k1,k2] += -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
        # a2_i_j+1
        M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
        # a2_i+1_j-1
        #M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        M[k1,lin(i+1,j,Nx)+Ntot] += -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j
        M[k1,lin(i+1,j,Nx)+Ntot] += Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j+1
        M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)

        # 2nd eq; rhs = 0
        # for a1 part
        # a1_i-1_j-1
        #M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        M[k2,lin(i-1,j,Nx)] += 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i-1_j
        M[k2,lin(i-1,j,Nx)] += Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
        # a1_i-1_j+1
        M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
        # a1_i_j-1
        #M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
        M[k2,lin(i,j,Nx)] += 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
        # a1_i,j
        M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
        # a1_i_j+1
        M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
        # a1_i+1_j-1
        #M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        M[k2,lin(i+1,j,Nx)] += -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i+1,j
        M[k2,lin(i+1,j,Nx)] += Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i+1_j+1
        M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
        
        # for a2 part
        # a2_i-1_j-1
        #M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        M[k2,lin(i-1,j,Nx)+Ntot] += 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i-1_j
        M[k2,lin(i-1,j,Nx)+Ntot] += Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
        # a2_i-1_j+1
        M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
        # a2_i_j-1
        #M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
        M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
        # a2_i_j
        M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
        # a2_i_j+1
        M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
        # a2_i+1_j-1
        #M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        M[k2,lin(i+1,j,Nx)+Ntot] += -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j
        M[k2,lin(i+1,j,Nx)+Ntot] += Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j+1
        M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)

    end # end j=1 copy, no corners

    # copy the j=Ny to j=Ny+1; no corners
    j=Ny
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
        P121_i_j_phalf = P121[i,j]
        Q121_i_j_mhalf = 0.5*(Q121[i,j-1] + Q121[i,j])
        Q121_i_j_phalf = Q121[i,j]
        P122_i_j_mhalf = 0.5*(P122[i,j-1] + P122[i,j])
        P122_i_j_phalf = P122[i,j]
        Q122_i_j_mhalf = 0.5*(Q122[i,j-1] + Q122[i,j])
        Q122_i_j_phalf = Q122[i,j]
        # sigma_22
        P221_i_j_mhalf = 0.5*(P221[i,j-1] + P221[i,j])
        P221_i_j_phalf = P221[i,j]
        Q221_i_j_mhalf = 0.5*(Q221[i,j-1] + Q221[i,j])
        Q221_i_j_phalf = Q221[i,j]
        P222_i_j_mhalf = 0.5*(P222[i,j-1] + P222[i,j])
        P222_i_j_phalf = P222[i,j]
        Q222_i_j_mhalf = 0.5*(Q222[i,j-1] + Q222[i,j])
        Q222_i_j_phalf = Q222[i,j]

        k = lin(i,j,Nx)
        k1 = k # for a1
        k2 = k + Ntot # for a2

        # 1st eq; rhs = 0
        # for a1 part
        # a1_i-1_j-1
        M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i-1_j
        M[k1,lin(i-1,j,Nx)] += Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
        # a1_i-1_j+1
        #M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
        M[k1,lin(i-1,j,Nx)] += -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
        # a1_i_j-1
        M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
        # a1_i,j
        M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
        # a1_i_j+1
        #M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
        M[k1,lin(i,j,Nx)] += 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
        # a1_i+1_j-1
        M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i+1,j
        M[k1,lin(i+1,j,Nx)] += Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
        # a1_i+1_j+1
        #M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
        M[k1,lin(i+1,j,Nx)] += 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
        
        # for a2 part
        # a2_i-1_j-1
        M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i-1_j
        M[k1,lin(i-1,j,Nx)+Ntot] += Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
        # a2_i-1_j+1
        #M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
        M[k1,lin(i-1,j,Nx)+Ntot] += -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
        # a2_i_j-1
        M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
        # a2_i_j
        M[k1,k2] += -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
        # a2_i_j+1
        #M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
        M[k1,lin(i,j,Nx)+Ntot] += 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
        # a2_i+1_j-1
        M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j
        M[k1,lin(i+1,j,Nx)+Ntot] += Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j+1
        #M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)
        M[k1,lin(i+1,j,Nx)+Ntot] += 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)

        # 2nd eq; rhs = 0
        # for a1 part
        # a1_i-1_j-1
        M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i-1_j
        M[k2,lin(i-1,j,Nx)] += Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
        # a1_i-1_j+1
        #M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
        M[k2,lin(i-1,j,Nx)] += -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
        # a1_i_j-1
        M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
        # a1_i,j
        M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
        # a1_i_j+1
        #M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
        M[k2,lin(i,j,Nx)] += 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
        # a1_i+1_j-1
        M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i+1,j
        M[k2,lin(i+1,j,Nx)] += Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
        # a1_i+1_j+1
        #M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
        M[k2,lin(i+1,j,Nx)] += 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
        
        # for a2 part
        # a2_i-1_j-1
        M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i-1_j
        M[k2,lin(i-1,j,Nx)+Ntot] += Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
        # a2_i-1_j+1
        #M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
        M[k2,lin(i-1,j,Nx)+Ntot] += -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
        # a2_i_j-1
        M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
        # a2_i_j
        M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
        # a2_i_j+1
        #M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
        M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
        # a2_i+1_j-1
        M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j
        M[k2,lin(i+1,j,Nx)+Ntot] += Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
        # a2_i+1_j+1
        #M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
        M[k2,lin(i+1,j,Nx)+Ntot] += 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)

    end # end copy the j=Ny to j=Ny+1; no corners


    #corner i=1,j=1
    j = 1
    i = 1
    # midpoints in x
    # sigma_11
    P111_i_mhalf_j = P111[i,j]
    P111_i_phalf_j = 0.5*(P111[i,j] + P111[i+1,j])
    Q111_i_mhalf_j = Q111[i,j]
    Q111_i_phalf_j = 0.5*(Q111[i,j] + Q111[i+1,j])
    P112_i_mhalf_j = P112[i,j]
    P112_i_phalf_j = 0.5*(P112[i,j] + P112[i+1,j])
    Q112_i_mhalf_j = Q112[i,j]
    Q112_i_phalf_j = 0.5*(Q112[i,j] + Q112[i+1,j])
    # sigma_12
    P121_i_mhalf_j = P121[i,j]
    P121_i_phalf_j = 0.5*(P121[i,j] + P121[i+1,j])
    Q121_i_mhalf_j = Q121[i,j]
    Q121_i_phalf_j = 0.5*(Q121[i,j] + Q121[i+1,j])
    P122_i_mhalf_j = P122[i,j]
    P122_i_phalf_j = 0.5*(P122[i,j] + P122[i+1,j])
    Q122_i_mhalf_j = Q122[i,j]
    Q122_i_phalf_j = 0.5*(Q122[i,j] + Q122[i+1,j])

    # midpoints in y
    # sigma_12
    P121_i_j_mhalf = P121[i,j]
    P121_i_j_phalf = 0.5*(P121[i,j] + P121[i,j+1])
    Q121_i_j_mhalf = Q121[i,j]
    Q121_i_j_phalf = 0.5*(Q121[i,j] + Q121[i,j+1])
    P122_i_j_mhalf = P122[i,j]
    P122_i_j_phalf = 0.5*(P122[i,j] + P122[i,j+1])
    Q122_i_j_mhalf = Q122[i,j]
    Q122_i_j_phalf = 0.5*(Q122[i,j] + Q122[i,j+1])
    # sigma_22
    P221_i_j_mhalf = P221[i,j]
    P221_i_j_phalf = 0.5*(P221[i,j] + P221[i,j+1])
    Q221_i_j_mhalf = Q221[i,j]
    Q221_i_j_phalf = 0.5*(Q221[i,j] + Q221[i,j+1])
    P222_i_j_mhalf = P222[i,j]
    P222_i_j_phalf = 0.5*(P222[i,j] + P222[i,j+1])
    Q222_i_j_mhalf = Q222[i,j]
    Q222_i_j_phalf = 0.5*(Q222[i,j] + Q222[i,j+1])

    k = lin(i,j,Nx)
    k1 = k # for a1
    k2 = k + Ntot # for a2

    # 1st eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_i0_j0 = f1_left[j]*(0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i-1_j
    #M[k1,lin(i-1,j,Nx)] = Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_i0_j1 = f1_left[j]*(Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy))
    # a1_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_i0_j2 = f1_left[j+1]*(-0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy))
    # a1_i_j-1
    #M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    M[k1,lin(i,j,Nx)] += 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)] += -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i+1,j
    M[k1,lin(i+1,j,Nx)] += Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i+1_j+1
    M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
    
    # for a2 part
    # a2_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_i0_j0 = f2_left[j]*(0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i-1_j
    #M[k1,lin(i-1,j,Nx)+Ntot] = Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_i0_j1 = f2_left[j]*(Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy))
    # a2_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_i0_j2 = f2_left[j+1]*(-0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy))
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

    # move to the rhs the i0 terms; pickup minus sign
    b[k1] = -eq1_a1_i0_j0 -eq1_a1_i0_j1 -eq1_a1_i0_j2 -eq1_a2_i0_j0 -eq1_a2_i0_j1 -eq1_a2_i0_j2

    # 2nd eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_i0_j0 = f1_left[j]*(0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i-1_j
    #M[k2,lin(i-1,j,Nx)] = Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_i0_j1 = f1_left[j]*(Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy))
    # a1_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_i0_j2 = f1_left[j+1]*(-0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy))
    # a1_i_j-1
    #M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    M[k2,lin(i,j,Nx)] += 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)] += -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i+1,j
    M[k2,lin(i+1,j,Nx)] += Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i+1_j+1
    M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
    
    # for a2 part
    # a2_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_i0_j0 = f2_left[j]*(0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i-1_j
    #M[k2,lin(i-1,j,Nx)+Ntot] = Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_i0_j1 = f2_left[j]*(Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy))
    # a2_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_i0_j2 = f2_left[j+1]*(-0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy))
    # a2_i_j-1
    #M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)+Ntot] += -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j
    M[k2,lin(i+1,j,Nx)+Ntot] += Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j+1
    M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
    
    # move to the rhs the i0 terms; pickup minus sign
    b[k2] = -eq2_a1_i0_j0 -eq2_a1_i0_j1 -eq2_a1_i0_j2 -eq2_a2_i0_j0 -eq2_a2_i0_j1 -eq2_a2_i0_j2

    # corner i=1, j=Ny
    i = 1 
    j = Ny
    # midpoints in x
    # sigma_11
    P111_i_mhalf_j = P111[i,j]
    P111_i_phalf_j = 0.5*(P111[i,j] + P111[i+1,j])
    Q111_i_mhalf_j = Q111[i,j]
    Q111_i_phalf_j = 0.5*(Q111[i,j] + Q111[i+1,j])
    P112_i_mhalf_j = P112[i,j]
    P112_i_phalf_j = 0.5*(P112[i,j] + P112[i+1,j])
    Q112_i_mhalf_j = Q112[i,j]
    Q112_i_phalf_j = 0.5*(Q112[i,j] + Q112[i+1,j])
    # sigma_12
    P121_i_mhalf_j = P121[i,j]
    P121_i_phalf_j = 0.5*(P121[i,j] + P121[i+1,j])
    Q121_i_mhalf_j = Q121[i,j]
    Q121_i_phalf_j = 0.5*(Q121[i,j] + Q121[i+1,j])
    P122_i_mhalf_j = P122[i,j]
    P122_i_phalf_j = 0.5*(P122[i,j] + P122[i+1,j])
    Q122_i_mhalf_j = Q122[i,j]
    Q122_i_phalf_j = 0.5*(Q122[i,j] + Q122[i+1,j])

    # midpoints in y
    # sigma_12
    P121_i_j_mhalf = 0.5*(P121[i,j-1] + P121[i,j])
    P121_i_j_phalf = P121[i,j]
    Q121_i_j_mhalf = 0.5*(Q121[i,j-1] + Q121[i,j])
    Q121_i_j_phalf = Q121[i,j]
    P122_i_j_mhalf = 0.5*(P122[i,j-1] + P122[i,j])
    P122_i_j_phalf = P122[i,j]
    Q122_i_j_mhalf = 0.5*(Q122[i,j-1] + Q122[i,j])
    Q122_i_j_phalf = Q122[i,j]
    # sigma_22
    P221_i_j_mhalf = 0.5*(P221[i,j-1] + P221[i,j])
    P221_i_j_phalf = P221[i,j]
    Q221_i_j_mhalf = 0.5*(Q221[i,j-1] + Q221[i,j])
    Q221_i_j_phalf = Q221[i,j]
    P222_i_j_mhalf = 0.5*(P222[i,j-1] + P222[i,j])
    P222_i_j_phalf = P222[i,j]
    Q222_i_j_mhalf = 0.5*(Q222[i,j-1] + Q222[i,j])
    Q222_i_j_phalf = Q222[i,j]

    k = lin(i,j,Nx)
    k1 = k # for a1
    k2 = k + Ntot # for a2

    # 1st eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_i0_jNym1 = f1_left[j-1]*(0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i-1_j
    #M[k1,lin(i-1,j,Nx)] = Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_i0_jNy = f1_left[j]*(Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy))
    # a1_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_i0_jNyp1 = f1_left[j]*(-0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy))
    # a1_i_j-1
    M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    #M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    M[k1,lin(i,j,Nx)] += 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i+1,j
    M[k1,lin(i+1,j,Nx)] += Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)] += 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
    
    # for a2 part
    # a2_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_i0_jNym1 = f2_left[j-1]*(0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i-1_j
    #M[k1,lin(i-1,j,Nx)+Ntot] = Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_i0_jNy = f2_left[j]*(Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy))
    # a2_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_i0_jNyp1 = f2_left[j]*(-0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy))
    # a2_i_j-1
    M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k1,k2] += -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    #M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
    M[k1,lin(i,j,Nx)+Ntot] += 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j
    M[k1,lin(i+1,j,Nx)+Ntot] += Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)+Ntot] += 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)

    # move to the rhs the i0 terms; pickup minus sign
    b[k1] = -eq1_a1_i0_jNym1 -eq1_a1_i0_jNy -eq1_a1_i0_jNyp1 -eq1_a2_i0_jNym1 -eq1_a2_i0_jNy -eq1_a2_i0_jNyp1

    # 2nd eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_i0_jNym1 = f1_left[j-1]*(0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i-1_j
    #M[k2,lin(i-1,j,Nx)] = Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_i0_jNy = f1_left[j]*(Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy))
    # a1_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_i0_jNyp1 = f1_left[j]*(-0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy))
    # a1_i_j-1
    M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    #M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    M[k2,lin(i,j,Nx)] += 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i+1,j
    M[k2,lin(i+1,j,Nx)] = Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)] += 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
    
    # for a2 part
    # a2_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_i0_jNym1 = f2_left[j-1]*(0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i-1_j
    #M[k2,lin(i-1,j,Nx)+Ntot] = Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_i0_jNy = f2_left[j]*(Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy))
    # a2_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_i0_jNyp1 = f2_left[j]*(-0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy))
    # a2_i_j-1
    M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    #M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j
    M[k2,lin(i+1,j,Nx)+Ntot] += Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)+Ntot] += 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)

    # move to the rhs the i0 terms; pickup minus sign
    b[k2] = -eq2_a1_i0_jNym1 -eq2_a1_i0_jNy -eq2_a1_i0_jNyp1 -eq2_a2_i0_jNym1 -eq2_a2_i0_jNy -eq2_a2_i0_jNyp1

    # corner i=Nx, j=1
    i = Nx
    j = 1
    # midpoints in x
    # sigma_11
    P111_i_mhalf_j = 0.5*(P111[i-1,j] + P111[i,j])
    P111_i_phalf_j = P111[i,j]
    Q111_i_mhalf_j = 0.5*(Q111[i-1,j] + Q111[i,j])
    Q111_i_phalf_j = Q111[i,j]
    P112_i_mhalf_j = 0.5*(P112[i-1,j] + P112[i,j])
    P112_i_phalf_j = P112[i,j]
    Q112_i_mhalf_j = 0.5*(Q112[i-1,j] + Q112[i,j])
    Q112_i_phalf_j = Q112[i,j]
    # sigma_12
    P121_i_mhalf_j = 0.5*(P121[i-1,j] + P121[i,j])
    P121_i_phalf_j = P121[i,j]
    Q121_i_mhalf_j = 0.5*(Q121[i-1,j] + Q121[i,j])
    Q121_i_phalf_j = Q121[i,j]
    P122_i_mhalf_j = 0.5*(P122[i-1,j] + P122[i,j])
    P122_i_phalf_j = P122[i,j]
    Q122_i_mhalf_j = 0.5*(Q122[i-1,j] + Q122[i,j])
    Q122_i_phalf_j = Q122[i,j]

    # midpoints in y
    # sigma_12
    P121_i_j_mhalf = P121[i,j]
    P121_i_j_phalf = 0.5*(P121[i,j] + P121[i,j+1])
    Q121_i_j_mhalf = Q121[i,j]
    Q121_i_j_phalf = 0.5*(Q121[i,j] + Q121[i,j+1])
    P122_i_j_mhalf = P122[i,j]
    P122_i_j_phalf = 0.5*(P122[i,j] + P122[i,j+1])
    Q122_i_j_mhalf = Q122[i,j]
    Q122_i_j_phalf = 0.5*(Q122[i,j] + Q122[i,j+1])
    # sigma_22
    P221_i_j_mhalf = P221[i,j]
    P221_i_j_phalf = 0.5*(P221[i,j] + P221[i,j+1])
    Q221_i_j_mhalf = Q221[i,j]
    Q221_i_j_phalf = 0.5*(Q221[i,j] + Q221[i,j+1])
    P222_i_j_mhalf = P222[i,j]
    P222_i_j_phalf = 0.5*(P222[i,j] + P222[i,j+1])
    Q222_i_j_mhalf = Q222[i,j]
    Q222_i_j_phalf = 0.5*(Q222[i,j] + Q222[i,j+1])

    k = lin(i,j,Nx)
    k1 = k # for a1
    k2 = k + Ntot # for a2

    # 1st eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)] += 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i-1_j
    M[k1,lin(i-1,j,Nx)] += Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
    # a1_i-1_j+1
    M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
    # a1_i_j-1
    #M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    M[k1,lin(i,j,Nx)] += 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_iNxp1_j0 = f1_right[j]*(-0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i+1,j
    #M[k1,lin(i+1,j,Nx)] = Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_iNxp1_j1 = f1_right[j]*(Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_iNxp1_j2 = f1_right[j+1]*(0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy))
    
    # for a2 part
    # a2_i-1_j-1
    #M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)+Ntot] += 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    # a2_i-1_j
    M[k1,lin(i-1,j,Nx)+Ntot] += Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
    # a2_i-1_j+1
    M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
    # a2_i_j-1
    #M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
    M[k1,lin(i,j,Nx)+Ntot] += 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k1,k2] += -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    #M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_iNxp1_j0 = f2_right[j]*(-0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j
    #M[k1,lin(i+1,j,Nx)+Ntot] = Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_iNxp1_j1 = f2_right[j]*(Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_iNxp1_j2 = f2_right[j+1]*(0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy))

    # move to the rhs the i0 terms; pickup minus sign
    b[k1] = -eq1_a1_iNxp1_j0 -eq1_a1_iNxp1_j1 -eq1_a1_iNxp1_j2 -eq1_a2_iNxp1_j0 -eq1_a2_iNxp1_j1 -eq1_a2_iNxp1_j2

    # 2nd eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)] += 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i-1_j
    M[k2,lin(i-1,j,Nx)] += Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
    # a1_i-1_j+1
    M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
    # a1_i_j-1
    #M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    M[k2,lin(i,j,Nx)] += 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_iNxp1_j0 = f1_right[j]*(-0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i+1,j
    #M[k2,lin(i+1,j,Nx)] = Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_iNxp1_j1 = f1_right[j]*(Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_iNxp1_j2 = f1_right[j+1]*(0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy))
    
    # for a2 part
    # a2_i-1_j-1
    #M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)+Ntot] += 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i-1_j
    M[k2,lin(i-1,j,Nx)+Ntot] += Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
    # a2_i-1_j+1
    M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
    # a2_i_j-1
    #M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_iNxp1_j0 = f2_right[j]*(-0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j
    #M[k2,lin(i+1,j,Nx)+Ntot] = Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_iNxp1_j1 = f2_right[j]*(Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_iNxp1_j2 = f2_right[j+1]*(0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy))
    
    # move to the rhs the i0 terms; pickup minus sign
    b[k2] = -eq2_a1_iNxp1_j0 -eq2_a1_iNxp1_j1 -eq2_a1_iNxp1_j2 -eq2_a2_iNxp1_j0 -eq2_a2_iNxp1_j1 -eq2_a2_iNxp1_j2

    # corner i=Nx, j=Ny
    i = Nx
    j = Ny
    # midpoints in x
    # sigma_11
    P111_i_mhalf_j = 0.5*(P111[i-1,j] + P111[i,j])
    P111_i_phalf_j = P111[i,j]
    Q111_i_mhalf_j = 0.5*(Q111[i-1,j] + Q111[i,j])
    Q111_i_phalf_j = Q111[i,j]
    P112_i_mhalf_j = 0.5*(P112[i-1,j] + P112[i,j])
    P112_i_phalf_j = P112[i,j]
    Q112_i_mhalf_j = 0.5*(Q112[i-1,j] + Q112[i,j])
    Q112_i_phalf_j = Q112[i,j]
    # sigma_12
    P121_i_mhalf_j = 0.5*(P121[i-1,j] + P121[i,j])
    P121_i_phalf_j = P121[i,j]
    Q121_i_mhalf_j = 0.5*(Q121[i-1,j] + Q121[i,j])
    Q121_i_phalf_j = Q121[i,j]
    P122_i_mhalf_j = 0.5*(P122[i-1,j] + P122[i,j])
    P122_i_phalf_j = P122[i,j]
    Q122_i_mhalf_j = 0.5*(Q122[i-1,j] + Q122[i,j])
    Q122_i_phalf_j = Q122[i,j]

    # midpoints in y
    # sigma_12
    P121_i_j_mhalf = 0.5*(P121[i,j-1] + P121[i,j])
    P121_i_j_phalf = P121[i,j]
    Q121_i_j_mhalf = 0.5*(Q121[i,j-1] + Q121[i,j])
    Q121_i_j_phalf = Q121[i,j]
    P122_i_j_mhalf = 0.5*(P122[i,j-1] + P122[i,j])
    P122_i_j_phalf = P122[i,j]
    Q122_i_j_mhalf = 0.5*(Q122[i,j-1] + Q122[i,j])
    Q122_i_j_phalf = Q122[i,j]
    # sigma_22
    P221_i_j_mhalf = 0.5*(P221[i,j-1] + P221[i,j])
    P221_i_j_phalf = P221[i,j]
    Q221_i_j_mhalf = 0.5*(Q221[i,j-1] + Q221[i,j])
    Q221_i_j_phalf = Q221[i,j]
    P222_i_j_mhalf = 0.5*(P222[i,j-1] + P222[i,j])
    P222_i_j_phalf = P222[i,j]
    Q222_i_j_mhalf = 0.5*(Q222[i,j-1] + Q222[i,j])
    Q222_i_j_phalf = Q222[i,j]

    k = lin(i,j,Nx)
    k1 = k # for a1
    k2 = k + Ntot # for a2

    # 1st eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    M[k1,lin(i-1,j-1,Nx)] = 0.25*(P111_i_mhalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    # a1_i-1_j
    M[k1,lin(i-1,j,Nx)] += Q111_i_mhalf_j/(Δx^2) + 0.25*(Q121_i_j_mhalf - Q121_i_j_phalf)/(Δx*Δy)
    # a1_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)] = -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)] += -0.25*(P111_i_mhalf_j + Q121_i_j_phalf)/(Δx*Δy)
    # a1_i_j-1
    M[k1,lin(i,j-1,Nx)] = 0.25*(P111_i_mhalf_j - P111_i_phalf_j)/(Δx*Δy) + P121_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k1,k1] += -(Q111_i_mhalf_j + Q111_i_phalf_j)/(Δx^2) -(P121_i_j_mhalf + P121_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    #M[k1,lin(i,j+1,Nx)] = 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    M[k1,lin(i,j,Nx)] += 0.25*(P111_i_phalf_j - P111_i_mhalf_j)/(Δx*Δy) + P121_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k1,lin(i+1,j-1,Nx)] = -0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_iNxp1_jNym1 = f1_right[j-1]*(-0.25*(P111_i_phalf_j + Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i+1,j
    #M[k1,lin(i+1,j,Nx)] = Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy)
    eq1_a1_iNxp1_jNy = f1_right[j]*(Q111_i_phalf_j/(Δx^2) + 0.25*(Q121_i_j_phalf - Q121_i_j_mhalf)/(Δx*Δy))
    # a1_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)] = 0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy)
    eq1_a1_iNxp1_jNyp1 = f1_right[j]*(0.25*(P111_i_phalf_j + Q121_i_j_phalf)/(Δx*Δy))
    
    # for a2 part
    # a2_i-1_j-1
    M[k1,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    # a2_i-1_j
    M[k1,lin(i-1,j,Nx)+Ntot] += Q112_i_mhalf_j/(Δx^2) + 0.25*(Q122_i_j_mhalf - Q122_i_j_phalf)/(Δx*Δy)
    # a2_i-1_j+1
    #M[k1,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)+Ntot] += -0.25*(P112_i_mhalf_j + Q122_i_j_phalf)/(Δx*Δy)
    # a2_i_j-1
    M[k1,lin(i,j-1,Nx)+Ntot] = 0.25*(P112_i_mhalf_j - P112_i_phalf_j)/(Δx*Δy) + P122_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k1,k2] += -(Q112_i_mhalf_j + Q112_i_phalf_j)/(Δx^2) -(P122_i_j_mhalf + P122_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    #M[k1,lin(i,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
    M[k1,lin(i,j,Nx)+Ntot] += 0.25*(P112_i_phalf_j - P112_i_mhalf_j)/(Δx*Δy) + P122_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    #M[k1,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_iNxp1_jNym1 = f2_right[j-1]*(-0.25*(P112_i_phalf_j + Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j
    #M[k1,lin(i+1,j,Nx)+Ntot] = Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy)
    eq1_a2_iNxp1_jNy = f2_right[j]*(Q112_i_phalf_j/(Δx^2) +0.25*(Q122_i_j_phalf - Q122_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j+1
    #M[k1,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy)
    eq1_a2_iNxp1_jNyp1 = f2_right[j]*(0.25*(P112_i_phalf_j + Q122_i_j_phalf)/(Δx*Δy))

    # move to the rhs the i0 terms; pickup minus sign
    b[k1] = -eq1_a1_iNxp1_jNym1 -eq1_a1_iNxp1_jNy -eq1_a1_iNxp1_jNyp1 -eq1_a2_iNxp1_jNym1 -eq1_a2_iNxp1_jNy -eq1_a2_iNxp1_jNyp1

    # 2nd eq; rhs = 0
    # for a1 part
    # a1_i-1_j-1
    M[k2,lin(i-1,j-1,Nx)] = 0.25*(P121_i_mhalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    # a1_i-1_j
    M[k2,lin(i-1,j,Nx)] += Q121_i_mhalf_j/(Δx^2) + 0.25*(Q221_i_j_mhalf - Q221_i_j_phalf)/(Δx*Δy)
    # a1_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)] = -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)] += -0.25*(P121_i_mhalf_j + Q221_i_j_phalf)/(Δx*Δy)
    # a1_i_j-1
    M[k2,lin(i,j-1,Nx)] = 0.25*(P121_i_mhalf_j - P121_i_phalf_j)/(Δx*Δy) + P221_i_j_mhalf/(Δy^2)
    # a1_i,j
    M[k2,k1] += -(Q121_i_mhalf_j + Q121_i_phalf_j)/(Δx^2) -(P221_i_j_mhalf + P221_i_j_phalf)/(Δy^2)
    # a1_i_j+1
    #M[k2,lin(i,j+1,Nx)] = 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    M[k2,lin(i,j,Nx)] += 0.25*(P121_i_phalf_j - P121_i_mhalf_j)/(Δx*Δy) + P221_i_j_phalf/(Δy^2)
    # a1_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)] = -0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_iNxp1_jNym1 = f1_right[j-1]*(-0.25*(P121_i_phalf_j + Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i+1,j
    #M[k2,lin(i+1,j,Nx)] = Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy)
    eq2_a1_iNxp1_jNy = f1_right[j]*(Q121_i_phalf_j/(Δx^2) + 0.25*(Q221_i_j_phalf - Q221_i_j_mhalf)/(Δx*Δy))
    # a1_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)] = 0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy)
    eq2_a1_iNxp1_jNyp1 = f1_right[j]*(0.25*(P121_i_phalf_j + Q221_i_j_phalf)/(Δx*Δy))
    
    # for a2 part
    # a2_i-1_j-1
    M[k2,lin(i-1,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    # a2_i-1_j
    M[k2,lin(i-1,j,Nx)+Ntot] += Q122_i_mhalf_j/(Δx^2) + 0.25*(Q222_i_j_mhalf - Q222_i_j_phalf)/(Δx*Δy)
    # a2_i-1_j+1
    #M[k2,lin(i-1,j+1,Nx)+Ntot] = -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)+Ntot] += -0.25*(P122_i_mhalf_j + Q222_i_j_phalf)/(Δx*Δy)
    # a2_i_j-1
    M[k2,lin(i,j-1,Nx)+Ntot] = 0.25*(P122_i_mhalf_j - P122_i_phalf_j)/(Δx*Δy) + P222_i_j_mhalf/(Δy^2)
    # a2_i_j
    M[k2,k2] += -(Q122_i_mhalf_j + Q122_i_phalf_j)/(Δx^2) -(P222_i_j_mhalf + P222_i_j_phalf)/(Δy^2)
    # a2_i_j+1
    #M[k2,lin(i,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    M[k2,lin(i,j,Nx)+Ntot] += 0.25*(P122_i_phalf_j - P122_i_mhalf_j)/(Δx*Δy) + P222_i_j_phalf/(Δy^2)
    # a2_i+1_j-1
    #M[k2,lin(i+1,j-1,Nx)+Ntot] = -0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_iNxp1_jNym1 = f2_right[j-1]*(-0.25*(P122_i_phalf_j + Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j
    #M[k2,lin(i+1,j,Nx)+Ntot] = Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy)
    eq2_a2_iNxp1_jNy = f2_right[j]*(Q122_i_phalf_j/(Δx^2) +0.25*(Q222_i_j_phalf - Q222_i_j_mhalf)/(Δx*Δy))
    # a2_i+1_j+1
    #M[k2,lin(i+1,j+1,Nx)+Ntot] = 0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy)
    eq2_a2_iNxp1_jNyp1 = f2_right[j]*(0.25*(P122_i_phalf_j + Q222_i_j_phalf)/(Δx*Δy))

    # move to the rhs the i0 terms; pickup minus sign
    b[k2] = -eq2_a1_iNxp1_jNym1 -eq2_a1_iNxp1_jNy -eq2_a1_iNxp1_jNyp1 -eq2_a2_iNxp1_jNym1 -eq2_a2_iNxp1_jNy -eq2_a2_iNxp1_jNyp1

    return M,b
end

function displacement_BC_left_right(t::Float64, p::Params)
    f1_left  = -1.0*t*1e-2*ones(p.Ny) # 0.0*ones(p.Ny)
    f2_left  = -0*1e-1*ones(p.Ny)
    
    f1_right =  0*t*1e-2*ones(p.Ny) # t*0.5*1e-1*ones(p.Ny)
    f2_right =  -0*1e-2*ones(p.Ny) # 0*1e-1*ones(p.Ny)
    return f1_left, f2_left, f1_right, f2_right
end