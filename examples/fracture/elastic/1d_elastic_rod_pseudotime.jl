using Serialization

# Define parameters
Base.@kwdef struct Params
    E::Float64 = 2.0        # Young's modulus
    gc::Float64 = 1e-2      # critical energy for fracture
    l::Float64 = 0.005      # regularization length (smaller to localize)
    k::Float64 = 1e-8       # small numerical factor
    L::Float64 = 1.0        # rod length
    N::Int = 2*512 + 1      # number of nodes (dx < l/3)
    dx::Float64 = L/(N-1)   # spatial step
    T::Int = 100 + 1        # number of steps for pseudotime
    dt::Float64 = 1/(T-1)   # step in pseudotime; when it is 1, we apply the full force on the right end
    tol::Float64 = 1e-6     # convergence tolerance
    max_iter::Int = 20      # max staggered iterations
    savefile::String = "fracture.jls" # file to save iterations
end

# Gaussian
function gaussian(x, x0, width, amp)
    return amp * exp.(-(x.-x0).^2 ./ (2*width^2))
end

# midpoint
function stagav(f::Vector)
    N = length(f)
    f_mid = zeros(N-1)
    for i in 1:N-1
        f_mid[i] = 0.5*(f[i]+f[i+1])
    end
    return f_mid
end

# strain u_x
function DuDx(u::Vector, p::Params)

    N = p.N
    idx = 1/p.dx
    u_x = zeros(N)

    # forward FD left
    u_x[1] = idx*(u[2]-u[1])

    # centred FD interior
    for i in 2:N-1
        u_x[i] = 0.5*idx*(u[i+1]-u[i-1])
    end

    # Neumann right
    #u_x[N] = F / ((E[N])*((1-d[N])^2 + p.k))

    # backward FD right
    u_x[N] = idx*(u[N]-u[N-1])
    return u_x
end

# solver for a single timestep
function ti_solver(u::Vector, d::Vector, H_old::Vector,
        ubar::Float64, E::Vector, p::Params)

    N = p.N
    l = p.l
    gc = p.gc
    idx = 1/p.dx
    u_x = DuDx(u,p)
    H = copy(H_old)

    for iter in 1:p.max_iter
        # copy to do convergence check at the end
        u_old = copy(u)
        d_old = copy(d)

        # Step 3: displacement
        g = (1.0 .- d).^2 .+ p.k
        g_mid = stagav(g)

        A = zeros(N,N)
        b = zeros(N)

        # Dirichlet left
        A[1,1] = 1.0
        b[1] = 0.0

        for i in 2:N-1
            if u_x[i] > 0
                A[i,i-1] = g_mid[i-1]*E[i]*idx
                A[i,i]   = -E[i]*idx*(g_mid[i]+g_mid[i-1])
                A[i,i+1] = g_mid[i]*E[i]*idx
            else
                A[i,i-1] = E[i]*idx
                A[i,i]   = -2*E[i]*idx
                A[i,i+1] = E[i]*idx
            end
        end

        # # Neumann right
        # A[N,N] = 1.0*idx
        # A[N,N-1] = -1.0*idx
        # b[N] = 2*F / (E[N]*g[N])
        """
        something is off with my Neuman BC; check again
        """

        # Dirichlet
        A[N,N] = 1.0
        b[N] = ubar

        u = A \ b
        
        # Step 1: history
        u_x = DuDx(u,p)
        H_new = 0.5*E.*max.(u_x,0).^2 # if u_x is negative we view it as compression and keep only tension
        H = max.(H,H_new) # keep the maximum at each spatial grid point
        """
        dont update the old history for every iteration of the same ti,
        only at the end of iter's where we have reached a good enough solution (convergence)
        """
        #H_old = copy(H) # update the old History

        println("   Iteration $iter: max(d_old) = ", maximum(d_old))

        # Step 2: phase-field
        M = zeros(N,N)
        f = 2*H

        for i in 2:N-1
            M[i,i-1] = -l*gc*idx^2
            M[i,i]   = gc/l + 2*l*gc*idx^2 + 2*H[i]
            M[i,i+1] = -l*gc*idx^2
        end

        # Neumann BC
        M[1,1] = gc/l + 2*l*gc*idx^2 + 2*H[1]#gc/l + 2*l*gc*idx^2 + 2*H[1]
        M[1,2] = -2*l*gc*idx^2#-2*l*gc*idx^2
        M[N,N] = gc/l + 2*l*gc*idx^2 + 2*H[N]#gc/l + 2*l*gc*idx^2 + 2*H[N]
        M[N,N-1] = -2*l*gc*idx^2#-2*l*gc*idx^2

        d = M \ f

        # convergence check
        res = maximum(abs.(u - u_old)) + maximum(abs.(d - d_old))
        if res < p.tol
            println("   res = $res | converged at iteration $iter")
            break
        end # end if

    end # end iters

    return u, d, H
end # end ti_solver

# apply ti_solver at each ti in the pseudotime domain
function pseudotime_solver(u0::Vector, d0::Vector, H0::Vector, 
    ubar::Float64, E::Vector, p::Params)

    N = p.N # number of spatial nodes
    T = p.T # number of pseudotemporal steps

    # pseudotime
    t = range(p.dt,1,T-1)
    # initiale matrices to save data
    u_all = zeros(N,T-1)
    d_all = zeros(N,T-1)
    H_all = zeros(N,T-1)

    # save ID
    d_all[:,1] = d0
    u_all[:,1] = u0
    H_all[:,1] = H0
    
    # prepare to first steps ti=1 that is t=0
    u = u0
    d = d0

    H_old = H0
    for ti in 2:T-1
        # prescribed displacement at ti
        ubar_i = t[ti]*ubar

        println("ti = $ti | ubar_i = $ubar_i:")
        # solve for ti
        u, d, H = ti_solver(u, d, H_old, ubar_i, E, p)
        H_old = copy(H)

        # save output at the end of iterations
        u_all[:,ti] = copy(u)
        d_all[:,ti] = copy(d)
        H_all[:,ti] = copy(H)
    end

    return u_all, d_all, H_all
end

function compute_H_from_d(d::Vector, H_old::Vector, p::Params)
    N = p.N
    l = p.l
    E = p.E
    gc = p.gc
    idx = 1/p.dx

    M = zeros(N,N)

    for i in 2:N-1
        M[i,i-1] = -l*gc*idx^2
        M[i,i]   = gc/l + 2*l*gc*idx^2 + 2*H_old[i]
        M[i,i+1] = -l*gc*idx^2
    end

    M[1,1] = gc/l + 2*l*gc*idx^2 + 2*H_old[1]
    M[1,2] = -2*l*gc*idx^2
    M[N,N] = gc/l + 2*l*gc*idx^2 + 2*H_old[N]
    M[N,N-1] = -2*l*gc*idx^2

    # Solve for H that satisfies the phase-field equation
    # 2H = M * d  =>  H = 0.5 * M * d
    H0 = 0.5 * abs.(M * d)

    return H0
end

# Run example
p = Params()
x = range(0,p.L,p.N)

u0 = zeros(p.N)
d0 = 0.01*rand(p.N)#gaussian(x, p.L/2, 0.05, 0.1)#zeros(p.N)
H0 = compute_H_from_d(d0, zeros(p.N), p)#0.1*gaussian(x, p.L/2, 0.05, 0.2)#0.1*rand(p.N)

# initial small defect
# defect_i = round(Int,p.N/2)
# d0[defect_i] = 0.1 
#u0 = gaussian(x, p.L/2, 0.5, 0.1)
# Define variable E with central soft spot possibly
E = ones(p.N) .* p.E
#E[defect_i-2:defect_i+2] .= p.E*0.5  # soften 5 points at center

# desired force at the right end
ubar = 1.2
u_sol, d_sol, H_sol = pseudotime_solver(u0, d0, H0, ubar, E, p)

# save all iterations
serialize(p.savefile, (x, u_sol, d_sol, H_sol))
