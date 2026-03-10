# # Build Q,P coefficients after linearization
# function build_QP(d,a1_bar,a2_bar,p::Params)

#     Nx = p.Nx
#     Ny = p.Ny
#     λ,μ = p.λ,p.μ

#     Q111 = zeros(Nx,Ny)
#     P111 = zeros(Nx,Ny)
#     Q112 = zeros(Nx,Ny)
#     P112 = zeros(Nx,Ny)

#     Q121 = zeros(Nx,Ny)
#     P121 = zeros(Nx,Ny)
#     Q122 = zeros(Nx,Ny)
#     P122 = zeros(Nx,Ny)

#     Q221 = zeros(Nx,Ny)
#     P221 = zeros(Nx,Ny)
#     Q222 = zeros(Nx,Ny)
#     P222 = zeros(Nx,Ny)

#     for j in 1:Ny
#         for i in 1:Nx

#             a1barx  = Dx(a1_bar,i,j,p)
#             a1bary  = Dy(a1_bar,i,j,p)
#             a2barx  = Dx(a2_bar,i,j,p)
#             a2bary  = Dy(a2_bar,i,j,p)

#             # # pressure
#             # pres = -(λ + 2μ/3)*(a1barx + a2bary)
#             # # mask
#             # if pres <= 0 
#             #     m_p = 1.0
#             # else    
#             #     m_p = g(d, i, j, p)
#             # end

#             # # from linearisation
#             # fbar = sqrt((a1barx - a2bary)^2 + (a1bary + a2barx)^2)

#             # if fbar == 0.0 # strictly
#             #     r1x = 1/6 - 1/(2^(3/2))
#             #     r1y = -1/(2^(3/2))
#             #     r2x = -1/(2^(3/2))
#             #     r2y = 1/6 + 1/(2^(3/2))
#             # else
#             #     r1x = 1/6 - (a1barx - a2bary)/(2fbar)
#             #     r1y = -(a1bary + a2barx)/(2fbar)
#             #     r2x = -(a1bary + a2barx)/(2fbar)
#             #     r2y = 1/6 + (a1barx - a2bary)/(2fbar)
#             # end

#             # q1 = (a1barx + a2bary)/6 - 0.5*fbar
#             # q2 = (a1barx + a2bary)/6 + 0.5*fbar
#             # # masks
#             # m_q1 = q1 <= 0 ? 1.0 : g(d, i, j, p)
#             # m_q2 = q2 <= 0 ? 1.0 : g(d, i, j, p)

#             # coefficients
#             # for a1x
#             Q111[i,j] = (λ + 2*μ)*g(d, i, j, p)
#             # for a1y
#             #P111[i,j] = 
#             # for a2x 
#             #Q112[i,j] = 
#             # for a2y
#             P112[i,j] = λ*g(d, i, j, p)

#             #Q121[i,j] = 
#             P121[i,j] = μ*g(d, i, j, p)
#             Q122[i,j] = μ*g(d, i, j, p)
#             #P122[i,j] = 

#             Q221[i,j] = λ*g(d, i, j, p)
#             #P221[i,j] = 
#             #Q222[i,j] = 
#             P222[i,j] = (λ + 2*μ)*g(d, i, j, p)
#         end
#     end

#     return Q111,P111,Q112, P112,
#            Q121,P121,Q122,P122,
#            Q221,P221,Q222,P222
# end


# Build Q,P coefficients after linearization
function build_QP(d,a1_bar,a2_bar,p::Params)

    Nx = p.Nx
    Ny = p.Ny
    λ,μ = p.λ,p.μ

    Q111 = zeros(Nx,Ny)
    P111 = zeros(Nx,Ny)
    Q112 = zeros(Nx,Ny)
    P112 = zeros(Nx,Ny)

    Q121 = zeros(Nx,Ny)
    P121 = zeros(Nx,Ny)
    Q122 = zeros(Nx,Ny)
    P122 = zeros(Nx,Ny)

    Q221 = zeros(Nx,Ny)
    P221 = zeros(Nx,Ny)
    Q222 = zeros(Nx,Ny)
    P222 = zeros(Nx,Ny)

    for j in 1:Ny
        for i in 1:Nx

            a1barx  = Dx(a1_bar,i,j,p)
            a1bary  = Dy(a1_bar,i,j,p)
            a2barx  = Dx(a2_bar,i,j,p)
            a2bary  = Dy(a2_bar,i,j,p)

            # pressure
            pres = (a1barx + a2bary)
            # mask
            if pres <= 0 
                mp = 1.0 + p.γstar*(1.0 - g(d, i, j, p))
            else    
                mp = g(d, i, j, p)
            end

            # # from linearisation
            # fbar = sqrt((a1barx - a2bary)^2 + (a1bary + a2barx)^2)

            # if fbar == 0.0 # strictly
            #     r1x = 1/6 - 1/(2^(3/2))
            #     r1y = -1/(2^(3/2))
            #     r2x = -1/(2^(3/2))
            #     r2y = 1/6 + 1/(2^(3/2))
            # else
            #     r1x = 1/6 - (a1barx - a2bary)/(2fbar)
            #     r1y = -(a1bary + a2barx)/(2fbar)
            #     r2x = -(a1bary + a2barx)/(2fbar)
            #     r2y = 1/6 + (a1barx - a2bary)/(2fbar)
            # end

            # q1 = (a1barx + a2bary)/6 - 0.5*fbar
            # q2 = (a1barx + a2bary)/6 + 0.5*fbar
            # # masks
            # m_q1 = q1 <= 0 ? 1.0 : g(d, i, j, p)
            # m_q2 = q2 <= 0 ? 1.0 : g(d, i, j, p)

            # coefficients
            # for a1x
            Q111[i,j] = (λ + 2*μ/3)*mp + 4*g(d,i,j,p)*μ/3 #λ + 2*μ
            # for a1y
            #P111[i,j] = 
            # for a2x 
            Q112[i,j] = (λ + 2*μ/3)*mp
            # for a2y
            P112[i,j] = -2*g(d,i,j,p)*μ/3 #λ

            #Q121[i,j] = 
            P121[i,j] = g(d,i,j,p)*μ #μ 
            Q122[i,j] = g(d,i,j,p)*μ #μ
            #P122[i,j] = 

            Q221[i,j] = (λ + 2*μ/3)*mp - 2*g(d,i,j,p)*μ/3 # λ
            #P221[i,j] = 
            Q222[i,j] = (λ + 2*μ/3)*mp
            P222[i,j] = 4*g(d,i,j,p)*μ/3 # λ + 2*μ
        end
    end

    return Q111,P111,Q112, P112,
           Q121,P121,Q122,P122,
           Q221,P221,Q222,P222
end
