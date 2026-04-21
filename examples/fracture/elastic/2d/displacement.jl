function build_displacement_system(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p::Params,
    f1_left, f1_right,
    f2_left, f2_right)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny
    Ndof = 2*Ntot

    M = zeros(Ndof,Ndof)
    b = zeros(Ndof)

    for j in 1:Ny
        for i in 2:Nx-1

    # build eq 1: 1/(Δx Δy)*(∫_E σ_11 Δy - ∫_W σ_11 Δy + ∫_N σ_12 Δx - ∫_S σ_12 Δx) = 0
    # add all the contributions to the matrix M and to the rhs b; mutates them
    eq1_int_E_σ11_Δy_over_volume!(
        Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq1_minus_int_W_σ11_Δy_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq1_int_N_σ12_Δx_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq1_minus_int_S_σ12_Δx_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    # build eq 2: 1/(Δx Δy)*(∫_E σ_12 Δy - ∫_W σ_12 Δy + ∫_N σ_22 Δx - ∫_S σ_22 Δx) = 0
    # add all the contributions to the matrix M and to the rhs b; mutates them
    eq2_int_E_σ12_Δy_over_volume!(
        Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq2_minus_int_W_σ12_Δy_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq2_int_N_σ22_Δx_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

    eq2_minus_int_S_σ22_Δx_over_volume!(Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222,
        p, M, b, i, j)

        end # end i
    end # end j

    # =========================
    # Dirichlet BC (left/right)
    # =========================

    # left boundary
    for j in 1:Ny
        k = lin(1,j,Nx)
        k1 = k
        k2 = k + Ntot
        M[k1,k1] = 1; b[k1] = f1_left[j]
        M[k2,k2] = 1; b[k2] = f2_left[j]
    end

    # right boundary
    for j in 1:Ny
        k = lin(Nx,j,Nx)
        k1 = k
        k2 = k + Ntot
        M[k1,k1] = 1; b[k1] = f1_right[j]
        M[k2,k2] = 1; b[k2] = f2_right[j]
    end

    return M,b
end

"""
1st eq contributions
"""
function eq1_int_E_σ11_Δy_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny

    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in x
    # =========================

    P111_i_phalf_j = 0.5*(P111[i,j]   + P111[i+1,j])
    Q111_i_phalf_j = 0.5*(Q111[i,j]   + Q111[i+1,j])
    P112_i_phalf_j = 0.5*(P112[i,j]   + P112[i+1,j])
    Q112_i_phalf_j = 0.5*(Q112[i,j]   + Q112[i+1,j])

    # add contributions
    
    # a1
    # x-derivative on east face
    M[k1,k1] += -Q111_i_phalf_j/(Δx^2)
    M[k1,lin(i+1,j,Nx)] += Q111_i_phalf_j/(Δx^2)
    # y derivative on east face
    M[k1,lin(i+1,jp,Nx)] += 0.25*P111_i_phalf_j/(Δx*Δy)
    M[k1,lin(i,jp,Nx)] += 0.25*P111_i_phalf_j/(Δx*Δy)
    M[k1,lin(i+1,jm,Nx)] += -0.25*P111_i_phalf_j/(Δx*Δy)
    M[k1,lin(i,jm,Nx)] += -0.25*P111_i_phalf_j/(Δx*Δy)

    # a2
    # x-derivative on east face
    M[k1,k2] += -Q112_i_phalf_j/(Δx^2)
    M[k1,lin(i+1,j,Nx)+Ntot] += Q112_i_phalf_j/(Δx^2)
    # y derivative on east face
    M[k1,lin(i+1,jp,Nx)+Ntot] += 0.25*P112_i_phalf_j/(Δx*Δy)
    M[k1,lin(i,jp,Nx)+Ntot] += 0.25*P112_i_phalf_j/(Δx*Δy)
    M[k1,lin(i+1,jm,Nx)+Ntot] += -0.25*P112_i_phalf_j/(Δx*Δy)
    M[k1,lin(i,jm,Nx)+Ntot] += -0.25*P112_i_phalf_j/(Δx*Δy)

end # end function

function eq1_minus_int_W_σ11_Δy_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny
    
    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in x
    # =========================

    P111_i_mhalf_j = 0.5*(P111[i-1,j] + P111[i,j])
    Q111_i_mhalf_j = 0.5*(Q111[i-1,j] + Q111[i,j])
    P112_i_mhalf_j = 0.5*(P112[i-1,j] + P112[i,j])
    Q112_i_mhalf_j = 0.5*(Q112[i-1,j] + Q112[i,j])

    # add contributions
    
    # a1
    # x-derivative on west face
    M[k1,k1] += -Q111_i_mhalf_j/(Δx^2)
    M[k1,lin(i-1,j,Nx)] += Q111_i_mhalf_j/(Δx^2)
    # y derivative on west face
    M[k1,lin(i-1,jp,Nx)] += -0.25*P111_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i,jp,Nx)] += -0.25*P111_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i-1,jm,Nx)] += 0.25*P111_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i,jm,Nx)] += 0.25*P111_i_mhalf_j/(Δx*Δy)

    # a2
    # x-derivative on west face
    M[k1,k2] += -Q112_i_mhalf_j/(Δx^2)
    M[k1,lin(i-1,j,Nx)+Ntot] += Q112_i_mhalf_j/(Δx^2)
    # y derivative on west face
    M[k1,lin(i-1,jp,Nx)+Ntot] += -0.25*P112_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i,jp,Nx)+Ntot] += -0.25*P112_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i-1,jm,Nx)+Ntot] += 0.25*P112_i_mhalf_j/(Δx*Δy)
    M[k1,lin(i,jm,Nx)+Ntot] += 0.25*P112_i_mhalf_j/(Δx*Δy)

end # end function

function eq1_int_N_σ12_Δx_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny
    
    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in y (with Neumann)
    # =========================

    P121_i_j_phalf = 0.5*(P121[i,j]  + P121[i,jp])
    Q121_i_j_phalf = 0.5*(Q121[i,j]  + Q121[i,jp])
    P122_i_j_phalf = 0.5*(P122[i,j]  + P122[i,jp])
    Q122_i_j_phalf = 0.5*(Q122[i,j]  + Q122[i,jp])

    # add contributions
    
    # a1
    # x-derivative on north face
    M[k1,k1] += -P121_i_j_phalf/(Δy^2)
    M[k1,lin(i,jp,Nx)] += P121_i_j_phalf/(Δy^2)
    # y derivative on north face
    M[k1,lin(i+1,jp,Nx)] += 0.25*Q121_i_j_phalf/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)] += 0.25*Q121_i_j_phalf/(Δx*Δy)
    M[k1,lin(i-1,jp,Nx)] += -0.25*Q121_i_j_phalf/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)] += -0.25*Q121_i_j_phalf/(Δx*Δy)

    # a2
    # x-derivative on north face
    M[k1,k2] += -P122_i_j_phalf/(Δy^2)
    M[k1,lin(i,jp,Nx)+Ntot] += P122_i_j_phalf/(Δy^2)
    # y derivative on north face
    M[k1,lin(i+1,jp,Nx)+Ntot] += 0.25*Q122_i_j_phalf/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)+Ntot] += 0.25*Q122_i_j_phalf/(Δx*Δy)
    M[k1,lin(i-1,jp,Nx)+Ntot] += -0.25*Q122_i_j_phalf/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)+Ntot] += -0.25*Q122_i_j_phalf/(Δx*Δy)

end # end function

function eq1_minus_int_S_σ12_Δx_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny

    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1

    # =========================
    # Midpoints in y (with Neumann)
    # =========================

    P121_i_j_mhalf = 0.5*(P121[i,jm] + P121[i,j])
    Q121_i_j_mhalf = 0.5*(Q121[i,jm] + Q121[i,j])
    P122_i_j_mhalf = 0.5*(P122[i,jm] + P122[i,j])
    Q122_i_j_mhalf = 0.5*(Q122[i,jm] + Q122[i,j])

    # add contributions
    
    # a1
    # x-derivative on south face
    M[k1,k1] += -P121_i_j_mhalf/(Δy^2)
    M[k1,lin(i,jm,Nx)] += P121_i_j_mhalf/(Δy^2)
    # y derivative on south face
    M[k1,lin(i+1,jm,Nx)] += -0.25*Q121_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)] += -0.25*Q121_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i-1,jm,Nx)] += 0.25*Q121_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)] += 0.25*Q121_i_j_mhalf/(Δx*Δy)

    # a2
    # x-derivative on south face
    M[k1,k2] += -P122_i_j_mhalf/(Δy^2)
    M[k1,lin(i,jm,Nx)+Ntot] += P122_i_j_mhalf/(Δy^2)
    # y derivative on south face
    M[k1,lin(i+1,jm,Nx)+Ntot] += -0.25*Q122_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i+1,j,Nx)+Ntot] += -0.25*Q122_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i-1,jm,Nx)+Ntot] += 0.25*Q122_i_j_mhalf/(Δx*Δy)
    M[k1,lin(i-1,j,Nx)+Ntot] += 0.25*Q122_i_j_mhalf/(Δx*Δy)

end # end function

"""
2nd eq contributions
"""
function eq2_int_E_σ12_Δy_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny

    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot


    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in x
    # =========================

    P121_i_phalf_j = 0.5*(P121[i,j]   + P121[i+1,j])
    Q121_i_phalf_j = 0.5*(Q121[i,j]   + Q121[i+1,j])
    P122_i_phalf_j = 0.5*(P122[i,j]   + P122[i+1,j])
    Q122_i_phalf_j = 0.5*(Q122[i,j]   + Q122[i+1,j])

    # add contributions
    
    # a1
    # x-derivative on east face
    M[k2,k1] += -Q121_i_phalf_j/(Δx^2)
    M[k2,lin(i+1,j,Nx)] += Q121_i_phalf_j/(Δx^2)
    # y-derivative on east face
    M[k2,lin(i+1,jp,Nx)] += 0.25*P121_i_phalf_j/(Δx*Δy)
    M[k2,lin(i,jp,Nx)] += 0.25*P121_i_phalf_j/(Δx*Δy)
    M[k2,lin(i+1,jm,Nx)] += -0.25*P121_i_phalf_j/(Δx*Δy)
    M[k2,lin(i,jm,Nx)] += -0.25*P121_i_phalf_j/(Δx*Δy)

    # a2
    # x-derivative on east face
    M[k2,k2] += -Q122_i_phalf_j/(Δx^2)
    M[k2,lin(i+1,j,Nx)+Ntot] += Q122_i_phalf_j/(Δx^2)
    # y derivative on east face
    M[k2,lin(i+1,jp,Nx)+Ntot] += 0.25*P122_i_phalf_j/(Δx*Δy)
    M[k2,lin(i,jp,Nx)+Ntot] += 0.25*P122_i_phalf_j/(Δx*Δy)
    M[k2,lin(i+1,jm,Nx)+Ntot] += -0.25*P122_i_phalf_j/(Δx*Δy)
    M[k2,lin(i,jm,Nx)+Ntot] += -0.25*P122_i_phalf_j/(Δx*Δy)

end # end function

function eq2_minus_int_W_σ12_Δy_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny

    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in x
    # =========================

    P121_i_mhalf_j = 0.5*(P121[i-1,j] + P121[i,j])
    Q121_i_mhalf_j = 0.5*(Q121[i-1,j] + Q121[i,j])
    P122_i_mhalf_j = 0.5*(P122[i-1,j] + P122[i,j])
    Q122_i_mhalf_j = 0.5*(Q122[i-1,j] + Q122[i,j])

    # add contributions
    
    # a1
    # x-derivative on west face
    M[k2,k1] += -Q121_i_mhalf_j/(Δx^2)
    M[k2,lin(i-1,j,Nx)] += Q121_i_mhalf_j/(Δx^2)
    # y derivative on west face
    M[k2,lin(i-1,jp,Nx)] += -0.25*P121_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i,jp,Nx)] += -0.25*P121_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i-1,jm,Nx)] += 0.25*P121_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i,jm,Nx)] += 0.25*P121_i_mhalf_j/(Δx*Δy)

    # a2
    # x-derivative on west face
    M[k2,k2] += -Q122_i_mhalf_j/(Δx^2)
    M[k2,lin(i-1,j,Nx)+Ntot] += Q122_i_mhalf_j/(Δx^2)
    # y derivative on west face
    M[k2,lin(i-1,jp,Nx)+Ntot] += -0.25*P122_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i,jp,Nx)+Ntot] += -0.25*P122_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i-1,jm,Nx)+Ntot] += 0.25*P122_i_mhalf_j/(Δx*Δy)
    M[k2,lin(i,jm,Nx)+Ntot] += 0.25*P122_i_mhalf_j/(Δx*Δy)

end # end function

function eq2_int_N_σ22_Δx_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny
    
    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jp = (j == Ny) ? j : j+1

    # =========================
    # Midpoints in y (with Neumann)
    # =========================

    P221_i_j_phalf = 0.5*(P221[i,j]  + P221[i,jp])
    Q221_i_j_phalf = 0.5*(Q221[i,j]  + Q221[i,jp])
    P222_i_j_phalf = 0.5*(P222[i,j]  + P222[i,jp])
    Q222_i_j_phalf = 0.5*(Q222[i,j]  + Q222[i,jp])

    # add contributions
    
    # a1
    # x-derivative on north face
    M[k2,k1] += -P221_i_j_phalf/(Δy^2)
    M[k2,lin(i,jp,Nx)] += P221_i_j_phalf/(Δy^2)
    # y derivative on north face
    M[k2,lin(i+1,jp,Nx)] += 0.25*Q221_i_j_phalf/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)] += 0.25*Q221_i_j_phalf/(Δx*Δy)
    M[k2,lin(i-1,jp,Nx)] += -0.25*Q221_i_j_phalf/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)] += -0.25*Q221_i_j_phalf/(Δx*Δy)

    # a2
    # x-derivative on north face
    M[k2,k2] += -P222_i_j_phalf/(Δy^2)
    M[k2,lin(i,jp,Nx)+Ntot] += P222_i_j_phalf/(Δy^2)
    # y derivative on north face
    M[k2,lin(i+1,jp,Nx)+Ntot] += 0.25*Q222_i_j_phalf/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)+Ntot] += 0.25*Q222_i_j_phalf/(Δx*Δy)
    M[k2,lin(i-1,jp,Nx)+Ntot] += -0.25*Q222_i_j_phalf/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)+Ntot] += -0.25*Q222_i_j_phalf/(Δx*Δy)

end # end function

function eq2_minus_int_S_σ22_Δx_over_volume!(
    Q111,P111,Q112,P112,
    Q121,P121,Q122,P122,
    Q221,P221,Q222,P222,
    p, M, b, i, j)

    Nx,Ny = p.Nx,p.Ny
    Δx,Δy = p.Δx,p.Δy
    Ntot = Nx*Ny

    k = lin(i,j,Nx)
    k1 = k
    k2 = k + Ntot

    # --- Neumann handling via mirroring ---
    jm = (j == 1)  ? j : j-1

    # =========================
    # Midpoints in y (with Neumann)
    # =========================

    P221_i_j_mhalf = 0.5*(P221[i,jm] + P221[i,j])
    Q221_i_j_mhalf = 0.5*(Q221[i,jm] + Q221[i,j])
    P222_i_j_mhalf = 0.5*(P222[i,jm] + P222[i,j])
    Q222_i_j_mhalf = 0.5*(Q222[i,jm] + Q222[i,j])

    # add contributions
    
    # a1
    # x-derivative on south face
    M[k2,k1] += -P221_i_j_mhalf/(Δy^2)
    M[k2,lin(i,jm,Nx)] += P221_i_j_mhalf/(Δy^2)
    # y derivative on south face
    M[k2,lin(i+1,jm,Nx)] += -0.25*Q221_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)] += -0.25*Q221_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i-1,jm,Nx)] += 0.25*Q221_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)] += 0.25*Q221_i_j_mhalf/(Δx*Δy)

    # a2
    # x-derivative on south face
    M[k2,k2] += -P222_i_j_mhalf/(Δy^2)
    M[k2,lin(i,jm,Nx)+Ntot] += P222_i_j_mhalf/(Δy^2)
    # y derivative on south face
    M[k2,lin(i+1,jm,Nx)+Ntot] += -0.25*Q222_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i+1,j,Nx)+Ntot] += -0.25*Q222_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i-1,jm,Nx)+Ntot] += 0.25*Q222_i_j_mhalf/(Δx*Δy)
    M[k2,lin(i-1,j,Nx)+Ntot] += 0.25*Q222_i_j_mhalf/(Δx*Δy)

end # end function