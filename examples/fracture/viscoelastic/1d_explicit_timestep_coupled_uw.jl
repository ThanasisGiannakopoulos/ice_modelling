using Serialization
using LinearAlgebra

# ============================================================
# Parameters
# ============================================================

Base.@kwdef struct Params
    # geometry
    L::Float64 = 1.0
    N::Int = 256
    dx::Float64 = L/(N-1)

    # time
    tf::Float64 = 1.0

    # phase field
    l::Float64 = 0.01
    kappa::Float64 = 1e-5

    # material
    ν::Float64 = 0.325
    A::Float64 = 1e-1
    n::Int = 3
    char_length::Float64 = 1e0
    char_displ::Float64 = 1e0
    Gc::Float64 = 1.0
    tau::Float64 = 1e-1*0.5

    # numerics
    tol::Float64 = 1e-9
    max_iter::Int = 20

    # critical energy for damage
    psi_crit::Float64 = 0.0#1e-2

    # for messages and save
    print_every::Int = 1
    save_every::Int = 1


    savefile::String = "run_explicit_adaptive_timestep.jls"
end

# ---------------------------
# Time dependent BC
# ---------------------------
function BC(t)
    return 2.0*t#0.1*abs(sin(t)))
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

# ============================================================
# For history field
# ============================================================

# get history from elastic displacement
function calculate_energy(u, w, p::Params)
    energy = zeros(length(u))
    ax = DfDx(u,p) .- DfDx(w,p)

    for i in eachindex(ax)
        # compute psi^+ depending on sign of ax
        if ax[i] > 0
            psi_plus = ax[i]^2 * (λ(p)/2 + 19*μ(p)/9)
        else
            psi_plus = ax[i]^2 * (8*μ(p)/9)
        end
        # use ⟨.⟩_+
        energy[i] = positive(psi_plus - p.psi_crit)
    end
    return energy
end

# get history from d; for initial data
function compute_H_from_d(d::Vector, p::Params)
    N = p.N
    l = p.l
    dx = p.dx
    C3 = CC3(p)#p.C3

    return abs.(d .- l^2 .*DfDxx(d,p))./(2*C3*l*(1 .- d .+ p.kappa))
end

# adaptive timestep using the nonlinear part
function get_dt(G, p)
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
# Coupled (u,d) solver for fixed w
# ============================================================

# you know w (fixed), get an u,d (fixed point iter), and History H
function get_udH(u_old, w_old, d_old, H_old, t, p::Params)
    N  = p.N
    dx = p.dx
    l = p.l 

    # BC
    T = BC(t)

    u = copy(u_old)
    d = copy(d_old)
    H = zeros(N)

    res = 0.0
    ii = 0
    for iter in 1:p.max_iter
        ii +=1
        # elastic strain
        a = u .- w_old
        ax = DfDx(a,p)

        # coefficients
        K = Kbar(ax,d,p)
        Kmid = stagav(K)

        # --------------------------------------------------
        # Allocate blocks
        # --------------------------------------------------
        A = zeros(N,N)
        a = zeros(N)

        # --------------------------------------------------
        # Boundary conditions (u)
        # --------------------------------------------------
        A[1,1] = 1.0
        a[1]   = 0.0

        A[N,N] = 1.0
        a[N] = T

        # --------------------------------------------------
        # Interior nodes
        # --------------------------------------------------
        for i in 2:N-1
            im = i-1
            ip = i+1

            # m=minus=i-1/2, p=plus=i+1/2
            km = Kmid[i-1]
            kp = Kmid[i]

            # ===== Momentum equation =====
            A[i,im] =  km#*idx^2
            A[i,i]  = -(km+kp)#*idx^2
            A[i,ip] =  kp#*idx^2

            a[i] = km*w_old[im]-(km+kp)*w_old[i] + kp*w_old[ip]

        end

        u_new = A \ a

        energy = calculate_energy(u_new,w_old,p)
        H = max.(H_old,energy)

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

        # Solve
        d_new = M \ f

        res = maximum(abs.(u_new-u)) + maximum(abs.(d_new-d))
        u .= u_new
        d .= d_new

        if res < p.tol
            break
        end
    end
    #println("   res = $(round(res, sigdigits=4)) | max(d) = $(round(maximum(d), digits=4)) | iter = $ii")

    return u, d, H, res, ii
end

# get f; dotwx=f, knowing u,w at the same timestep
function get_f(u,w,d,p::Params)

    # get f, dotwx=f
    a = u .- w
    ax = DfDx(a,p)
    F = Fcoef(ax,d,p)
    f = F.*ax.^p.n
    return f
end

# knowing the rhs, invert to get w
function get_w(w_old, rhs, dt, p::Params)
    N = p.N
    dx = p.dx
    w_new = zeros(N)
    
    # initiate matrix, vector for inversion w_new = A \ a
    A = zeros(N,N)
    a = zeros(N)

    # BC left, Dirichlet
    A[1,1] = 1.0
    a[1] = 0.0
    # BC right, Neumann, backward 1st order FD
    A[N,N] = 1.0
    A[N,N-1] = -1.0
    a[N] = w_old[N] - w_old[N-1] + 2.0*dt*dx*rhs[N]

    for i in 2:N-1
        A[i,i-1] = -1.0
        A[i,i+1] = 1.0
        a[i] = w_old[i+1] - w_old[i-1] + 2.0*dt*dx*rhs[i]
    end

    w_new = A \ a
    return w_new
end

# RK4 timestep, returns u,w,d,H at the next timestep, can be applied for adaptive dt
function RK4_timestep_uwdH(u_old, w_old, d_old, H_old, t, dt, p::Params)
    
    N = p.N

    # 1st step
    #println("get k1")
    k1 = zeros(N)
    k1 = get_f(u_old,w_old,d_old,p)

    # 2nd step
    #println("get k2")
    u1 = zeros(N)
    d1 = zeros(N)
    H1 = zeros(N) # redundant?

    w1 = w_old .+ k1*0.5*dt
    t1 = t + 0.5*dt
    u1, d1, H1, _ = get_udH(u_old, w1, d_old, H_old, t1, p)
    
    k2 = zeros(N)
    k2 = get_f(u1,w1,d1,p)

    #3rd step
    #println("get k3")
    u2 = zeros(N)
    d2 = zeros(N)
    H2 = zeros(N) # redundant?

    w2 = w_old .+ k2*0.5*dt
    t2 = t1
    # find again u,d,H for w2
    # not sure if it makes sense to use u1,d1,H1 instead
    u2, d2, H2, _ = get_udH(u_old, w2, d_old, H_old, t2, p)
    
    k3 = zeros(N)
    k3 = get_f(u2,w2,d2,p)

    # 4th step
    #println("get k4")
    u3 = zeros(N)
    d3 = zeros(N)
    H3 = zeros(N) # redundant?

    w3 = w_old .+ k3*dt
    t3 = t + dt
    # find again u,d,H for w2
    # not sure if it makes sense to use u1,d1,H1 instead
    u3, d3, H3, _ = get_udH(u_old, w3, d_old, H_old, t3, p)
    
    k4 = zeros(N)
    k4 = get_f(u3,w3,d3,p)

    # get w at the next timestep; rhs is an average i.e. rhs=(k1+2k2+2k3+k4)/6 )
    #println("get new sols")
    rhs = (k1 .+ 2*k2 .+ 2*k3 .+k4)/6
    w_new = get_w(w_old, rhs, dt, p)
    # now get u_new, d_new, H_new, for known w_new
    u_new, d_new, H_new, res, ii = get_udH(u_old, w_new, d_old, H_old, t+dt, p)

    return u_new, w_new, d_new, H_new, res, ii
end

# ============================================================
# Main time loop
# ============================================================

function run_simulation(p::Params)
    x = range(0,p.L,p.N)

    println("ν = $(p.ν) | μ = $(μ(p)) |  C2 = $(CC2(p)) | C3 = $(CC3(p))")

    u = zeros(p.N)
    w = zeros(p.N)
    d = gaussian(x, p.L/2, 0.02, 0.3)#0.01*rand(p.N)
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

        # get timestep with nonlinearities for current time
        a = u .- w
        ax = DfDx(a,p)
        G = Gcoef(ax, p)
        dt = get_dt(G,p)

        # calculate new sols
        u_new, w_new, d_new, H_new, res, ii = RK4_timestep_uwdH(u, w, d, H, t, dt, p)
        
        # roll-over new sols to active vars
        u = copy(u_new)
        w = copy(w_new)
        d = copy(d_new)
        H = copy(H_new)

        # update time
        t += dt
        step +=1

        # if res >p.tol
        #     println("Time step $step | t = $(round(t, digits=4))")
        #     append!(u_list, [copy(u)])
        #     append!(w_list, [copy(w)])
        #     append!(d_list, [copy(d)])
        #     append!(H_list, [copy(H)])
        #     append!(t_list, [copy(t)])
        #     serialize(p.savefile,(x,u_list,w_list,d_list,H_list, t_list))
        #     break
        # end
        if step % p.print_every == 0
            println("timestep = $step | t = $(round(t, digits=4)) | res = $(round(res, sigdigits=4)) | max(d) = $(round(maximum(d), digits=4)) | iter = $ii")
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
