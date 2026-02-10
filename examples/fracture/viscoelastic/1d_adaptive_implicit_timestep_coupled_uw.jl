using Serialization
using LinearAlgebra

# ============================================================
# Parameters
# ============================================================

Base.@kwdef struct Params
    # geometry
    L::Float64 = 10.0
    N::Int = 512
    dx::Float64 = L/(N-1)

    # time
    #NT::Int = 2*290#808
    #dt::Float64 = 1e-1*0.5*0.5*0.5
    tf::Float64 = 10.0

    # phase field
    l::Float64 = 0.1
    kappa::Float64 = 1e-5

    # material
    ν::Float64 = 0.325
    A::Float64 = 1e-2 #1.2*1e-25
    n::Int = 3
    #C2::Float64 = 11.82#0.25
    #C3::Float64 = 3520.75
    char_length::Float64 = 1e-1
    char_displ::Float64 = 1e-1
    Gc::Float64 = 0.6
    tau::Float64 = 1e-1*4

    # numerics
    tol::Float64 = 1e-9
    max_iter::Int = 20

    # critical energy for damage
    psi_crit::Float64 = 0.01#1e-2

    # for messages and save
    print_every::Int = 1
    save_every::Int = 1

    savefile::String = "run_adaptive_implicit_coupled_uw.jls"
end

# ============================================================
# Material properties
# ============================================================

# Young's modules tuned by Poisson ratio ν to satisfy mass balance
E(p) = 2*(1+p.ν)^2/(1-p.ν)

# Lame Parameters
λ(p) = E(p)*p.ν / ((1+p.ν)*(1-2p.ν))
μ(p) = E(p)/(2*(1+p.ν))

# model's dimless constants
CC2(p) = p.A^(1/p.n)*(p.char_displ/p.char_length)^(1-1/p.n)*μ(p)*p.tau^(1/p.n)
CC3(p) = μ(p)*p.char_displ^2/(p.Gc*p.char_length)

# ============================================================
# Utilities
# ============================================================

# Gaussian
function gaussian(x, x0, width, amp)
    return amp * exp.(-(x.-x0).^2 ./ (2*width^2))
end

# ⟨x⟩_+
positive(x) = 0.25 * (x + abs(x))

# centred FD with forward/backward on ends
function DfDx(f::Vector, p::Params)
    N = p.N
    dx = p.dx
    fx = zeros(N)
    fx[1] = (f[2]-f[1])/dx
    for i in 2:N-1
        fx[i] = (f[i+1]-f[i-1])/(2dx)
    end
    fx[N] = (f[N]-f[N-1])/dx
    return fx
end

# centred FD with forward/backward on ends
function DfDxx(f::Vector, p::Params)
    N = p.N
    dx = p.dx
    idx = 1/dx
    fxx = zeros(N)
    fxx[1] = (2*f[1]-5*f[2]+4*f[3]-f[4])*idx^2
    for i in 2:N-1
        fxx[i] = (f[i+1]-2*f[i]+f[i-1])*idx^2
    end
    fxx[N] = (2*f[N]-5*f[N-1]+4*f[N-2]-f[N-3])*idx^2
    return fxx
end

# midpoint
# stagav(f)_[i]=f_{i+1/2}, stagav(f)_[i-1]=f_{i-1/2} 
stagav(f) = 0.5 .* (f[1:end-1] .+ f[2:end])

# degradation function
g(d,p) = (1 .- d).^2 .+ p.kappa

# ============================================================
# Effective stiffness K̄(d, a_x), 
# with a_x = u_x - w_x (elastic displacement)
# ============================================================

function Kbar(ax, d, p)
    K0 = λ(p)/2 + μ(p)/3 + 4/3
    gval = g(d,p)
    K = similar(ax)
    for i in eachindex(ax)
        if ax[i] > 0
            K[i] = 2*gval[i]*K0
        else
            K[i] = 2*K0
        end
    end
    return K
end

# ============================================================
# Glen viscosity coefficients
# ============================================================

function Fcoef(ax, d, p)
    gval = g(d,p)
    n = p.n
    A = p.A
    C2 = CC2(p) #p.C2
    F = similar(d)
    for i in eachindex(d)
        if ax[i] > 0
            F[i] = 2^((1+n)/2)*A*
                   (C2*(16/9 + 8/(9*gval[i])))^n
        else
            F[i] = 2^((1+n)/2)*A*
                   (C2*(8/9 + 16/(9*gval[i])))^n
        end
    end
    return F
end

function Gcoef(ax, p)
    n = p.n
    return ax.^(n-1)
end

# adaptive timestep using the nonlinear part
function find_dt(G, p)
    dx = p.dx
    maxG = maximum(abs.(G))
    if maxG>1
        dt = 0.25*dx/maxG
    else
        dt = 0.25*dx
    end
    return dt
end

# ============================================================
# Coupled (u,w) solver for one timestep
# ============================================================

function solve_uw_step(u_old, w_old, d, T, p::Params)
    N  = p.N
    dx = p.dx
    idx = 1/dx
    n = p.n

    u = copy(u_old)
    w = copy(w_old)

    res = 0.0
    ii = 0

    # initiate dt s.t. it is changed and exported later
    dt = 0.25*dx

    for iter in 1:p.max_iter
        ii +=1
        # elastic strain
        a = u .- w
        ax = DfDx(a,p)

        # coefficients
        K = Kbar(ax,d,p)
        Kmid = stagav(K)

        F = Fcoef(ax,d,p)
        G = Gcoef(ax,p)
        Fmid = stagav(F)
        Gmid = stagav(G)
        Fx = DfDx(F, p)
        #Gx = DfDx(G, p)
        Gx = (n-1)*ax.^(n-2) .*DfDxx(a,p)

        # get timestep in each fixed point iteration
        dt = copy(find_dt(G,p))
        # println("dt = $(dt)")

        # --------------------------------------------------
        # Allocate blocks
        # --------------------------------------------------
        A = zeros(N,N)
        B = zeros(N,N)
        C = zeros(N,N)
        D = zeros(N,N)

        a = zeros(N)
        c = zeros(N)

        # --------------------------------------------------
        # Boundary conditions (u)
        # --------------------------------------------------
        A[1,1] = 1.0
        a[1]   = 0.0

        # # to prescribe traction T
        # A[N,N] = 1.0*idx
        # A[N,N-1] = -1.0*idx
        # B[N,N] = -1.0*idx
        # B[N,N-1] = 1.0*idx
        # a[N]   = T/K[N]

        A[N,N] = 1.0
        a[N] = T

        # --------------------------------------------------
        # Boundary conditions (w)
        # --------------------------------------------------
        # Dirichlet on left end; acts as gauge fixing since we evolve w_x
        D[1,1] = 1.0
        c[1] = 0.0 #0.0

        # Neuman with backward FD on right end
        D[N,N] = 1.0*idx
        D[N,N-1] = -1.0*idx
        ux_minus_wx_old = (u_old[N]-u_old[N-1])/dx - (w_old[N]-w_old[N-1])/dx
        c[N] = (w_old[N] - w_old[N-1])/dx + dt*F[N]*(ux_minus_wx_old)^3
        #c[N] = (w_old[N] - w_old[N-1])/dx + dt*F[N]*(T/K[N])^3#dt*F[N]*G[N]*(T/K[N])

        # --------------------------------------------------
        # Interior nodes
        # --------------------------------------------------
        for i in 2:N-1
            im = i-1
            ip = i+1

            # m=minus=i-1/2, p=plus=i+1/2
            km = Kmid[i-1]
            kp = Kmid[i]
            Fm = Fmid[i-1]
            Fp = Fmid[i]
            Gm = Gmid[i-1]
            Gp = Gmid[i]

            # ===== Momentum equation =====
            A[i,im] =  km#*idx^2
            A[i,i]  = -(km+kp)#*idx^2
            A[i,ip] =  kp#*idx^2

            B[i,im] = -km#*idx^2
            B[i,i]  =  (km+kp)#*idx^2
            B[i,ip] = -kp#*idx^2

            # ===== Viscosity equation =====
            C[i,im] =  0.5*dt*Fm*Gm
            C[i,i]  = -dt*(0.5*Fp*Gp-0.5*Fm*Gm-dx*(Fx[i]*G[i]+F[i]*Gx[i]))
            C[i,ip] = -0.5*dt*Fp*Gp

            D[i,im] = -0.5*(1 + dt*Fm*Gm)
            D[i,i]  = -dt*(0.5*Fp*Gp-0.5*Fm*Gm-dx*(Fx[i]*G[i]+F[i]*Gx[i]))
            D[i,ip] =  0.5*(1 + dt*Fp*Gp)

            # RHS
            c[i] = (w_old[ip] - w_old[im]) / 2
        end

        # --------------------------------------------------
        # Assemble full system
        # --------------------------------------------------
        M = [A  B;
             C  D]

        rhs = [a; c]

        sol = M \ rhs

        u_new = sol[1:N]
        w_new = sol[N+1:end]

        res = maximum(abs.(u_new-u)) + maximum(abs.(w_new-w))
        #println("   Iteration $iter: max(d) = ", maximum(d))
        u .= u_new
        w .= w_new

        if res < p.tol
            break
        end
    end
    # println("   res = $res | max(d) =  $(maximum(d)) | iter =  $iter")
    println("   res = $(round(res, sigdigits=4)) | max(d) = $(round(maximum(d), digits=4)) | iter = $ii | dt = $dt")


    return u, w, dt, res
end


# ============================================================
# Phase-field history
# ============================================================

# get history from elastic displacement
function update_history!(H, u, w, p::Params)
    ax = DfDx(u,p) .- DfDx(w,p)

    for i in eachindex(ax)
        # compute psi^+ depending on sign of ax
        if ax[i] > 0
            psi_plus = ax[i]^2 * (λ(p)/2 + 19*μ(p)/9)
        else
            psi_plus = ax[i]^2 * (8*μ(p)/9)
        end
        # use ⟨.⟩_+
        H[i] = max(H[i], positive(psi_plus - p.psi_crit))
    end
end

# get history from d; for initial data
function compute_H_from_d(d::Vector, p::Params)
    N = p.N
    l = p.l
    dx = p.dx
    C3 = CC3(p)#p.C3

    return abs.(d .- l^2 .*DfDxx(d,p))./(2*C3*l*(1 .- d .+ p.kappa))
end
# ============================================================
# Phase-field solver
# ============================================================

function solve_phasefield(H, p::Params)
    N = p.N
    dx = p.dx
    l = p.l
    C3 = CC3(p) #p.C3

    M = zeros(N,N)
    f = zeros(N)

    # Interior points
    for i in 2:N-1
        M[i, i-1] = -l^2 / dx^2
        M[i, i]   = 1 + 2*l^2/dx^2 + 2*C3*l*H[i]
        M[i, i+1] = -l^2 / dx^2
        f[i]      = 2*C3*l*H[i]
    end

    # Neumann BC at i=1 (ghost node)
    M[1,1] = 1 + 2*l^2/dx^2 + 2*C3*l*H[1]
    M[1,2] = -2*l^2/dx^2
    f[1]   = 2*C3*l*H[1]

    # Neumann BC at i=N (ghost node)
    M[N,N]   = 1 + 2*l^2/dx^2 + 2*C3*l*H[N]
    M[N,N-1] = -2*l^2/dx^2
    f[N]     = 2*C3*l*H[N]

    # M[1,1] = -1.0/dx
    # M[1,2] = 1/0/dx
    # M[N,N] = 1.0/dx
    # M[N,N-1] = -1.0/dx
    # Solve
    d_new = M \ f

    return d_new
end


# ---------------------------
# Time dependent BC
# ---------------------------
function traction(t, p::Params)
    return 2.0*t#0.1*abs(sin(t)))
end


# ============================================================
# Main time loop
# ============================================================

function run_simulation(p::Params)
    x = range(0,p.L,p.N)

    println("ν = $(p.ν) | μ = $(μ(p)) |  C2 = $(CC2(p)) | C3 = $(CC3(p))")

    u = zeros(p.N)
    w = zeros(p.N)
    d = gaussian(x, p.L/2, 0.2, 0.3)#0.01*rand(p.N)
    H = zeros(p.N)
    # initial history from initial damage
    H = compute_H_from_d(d, p)

    u_list = []
    w_list = []
    d_list = []
    H_list = []
    t_list = []

    t = 0.0
    step = 0
    
    append!(u_list, [copy(u)])
    append!(w_list, [copy(w)])
    append!(d_list, [copy(d)])
    append!(H_list, [copy(H)])
    append!(t_list, [copy(t)])

    println("Time step $step | t = $(round(t, digits=4))")
    println("max(d0) = $(maximum(d))")
    
    #for k in 2:p.NT
    while t < p.tf
        # BC for total displacement
        T = traction(t, p)

        u, w, dt, res = solve_uw_step(u,w,d,T,p)
        update_history!(H,u,w,p)
        d = solve_phasefield(H,p)

        t += dt
        step +=1

        if res >p.tol
            println("Time step $step | t = $(round(t, digits=4))")
            append!(u_list, [copy(u)])
            append!(w_list, [copy(w)])
            append!(d_list, [copy(d)])
            append!(H_list, [copy(H)])
            append!(t_list, [copy(t)])
            serialize(p.savefile,(x,u_list,w_list,d_list,H_list, t_list))
            break
        end
        if step % p.print_every == 0
            println("Time step $step | t = $(round(t, digits=4))")
        end

        if step % p.save_every == 0
            append!(u_list, [copy(u)])
            append!(w_list, [copy(w)])
            append!(d_list, [copy(d)])
            append!(H_list, [copy(H)])
            append!(t_list, [copy(t)])
        end
    end

    # save last timestep
    append!(u_list, [copy(u)])
    append!(w_list, [copy(w)])
    append!(d_list, [copy(d)])
    append!(H_list, [copy(H)])
    append!(t_list, [copy(t)])

    serialize(p.savefile,(x,u_list,w_list,d_list,H_list, t_list))
end

# ============================================================
# Run
# ============================================================

p = Params()
run_simulation(p)
