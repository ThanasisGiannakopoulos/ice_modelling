using Serialization

# parameters
Base.@kwdef struct Params
    E::Float64 = 2.0        # Young's modulus
    gc::Float64 = 0.2       # critical energy for fracture
    l::Float64 = 0.01       # regularization length
    k::Float64 = 1e-6       # small numerical factor
    L::Float64 = 1.0        # rod length
    N::Int = 512 + 1        # number of nodes
    dx::Float64 = L/(N-1)   # grid step
    tol::Float64 = 1e-6     # for convergence tolerance
    max_iter::Int = 20    # max iterations
    savefile::String = "fracture.jls" # file to save iterations
end

# midpoint
function stagav(f)
    N = length(f)
    f_mid = zeros(N-1)
    for i in 1:N-1
        f_mid[i] = 0.5*(f[i]+f[i+1])
    end
    return f_mid
end


# Gaussian
function gaussian(x, x0, width, amp)
    return amp * exp.(-(x.-x0).^2 ./ (2*width^2))
end
# Gaussian derivative; positive left, negative right of x0
function gaussian_derivative(x, x0, width, amp)
    return -amp * (x .- x0) ./ width^2 .* exp.(-(x .- x0).^2 ./ (2*width^2))
end

# compute u_x (with variable E)
function DuDx(u, d, F, E_vec, p::Params)
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
    u_x[N] = F / ((E_vec[N])*((1-d[N])^2 + p.k))

    return u_x
end

# main solver
function main_solver(x, u, d, H, F, p::Params)
    N = p.N
    idx = 1/p.dx
    H_old = H

    # Save arrays per iteration
    u_all = []
    d_all = []
    H_all = []

    # Define variable E
    E_vec = ones(N) .* p.E
    mid = round(Int, N/2)
    E_vec[mid-2:mid+2] .= p.E*0.5  # soften 5 points at center

    # initial save
    push!(u_all, copy(u))
    push!(d_all, copy(d))
    push!(H_all, copy(H_old))

    for iter = 1:p.max_iter
        u_old = copy(u)
        d_old = copy(d)

        # Step 1: history
        u_x = DuDx(u,d,F,E_vec,p)
        H_new = 0.5*E_vec.*max.(u_x,0).^2
        H = max.(H_old,H_new)
        H_old = copy(H)

        println("Iteration $iter: max(d) = ", maximum(d))

        # Step 2: phase-field
        M = zeros(N,N)
        f = 2*p.l*H

        for i in 2:N-1
            M[i,i-1] = -p.l*p.gc*idx^2
            M[i,i]   = p.gc/p.l + 2*p.l*p.gc*idx^2 + 2*H[i]
            M[i,i+1] = -p.l*p.gc*idx^2
        end

        # BC
        M[1,1] = p.gc/p.l + 2*p.l*p.gc*idx^2 + 2*H[1]
        M[1,2] = -2*p.l*p.gc*idx^2
        M[N,N] = p.gc/p.l + 2*p.l*p.gc*idx^2 + 2*H[N]
        M[N,N-1] = -2*p.l*p.gc*idx^2

        d = M \ f

        # Step 3: displacement
        g = (1.0 .- d).^2 .+ p.k
        g_mid = stagav(g)

        A = zeros(N,N)
        #b = zeros(N)
        b = -gaussian_derivative(x, 0.5, 0.05, 0.01)


        # Dirichlet left
        A[1,1] = 1.0
        b[1] = 0.0

        for i in 2:N-1
            A[i,i-1] = g_mid[i-1]*E_vec[i]*idx^2
            A[i,i]   = -E_vec[i]*idx^2*(g_mid[i]+g_mid[i-1])
            A[i,i+1] = g_mid[i]*E_vec[i]*idx^2
        end

        # Neumann right
        A[N,N] = 1.0
        A[N,N-1] = -1.0
        b[N] = 2*F / (E_vec[N]*g[N]) * p.dx

        u = A \ b

        # save iteration
        push!(u_all, copy(u))
        push!(d_all, copy(d))
        push!(H_all, copy(H))

        # convergence check
        res = maximum(abs.(u - u_old)) + maximum(abs.(d - d_old))
        if res < p.tol
            println("Converged at iteration $iter")
            break
        end
    end

    return u_all, d_all, H_all, E_vec
end

# Run example
p = Params()
x = range(0,p.L,p.N)

u = zeros(p.N)
d = zeros(p.N)
H = zeros(p.N)
# initial small defect
defect_i = round(Int,p.N/2)
d[defect_i] = 0.1 
# initial non-zero history where the defect is, st. it is not smoothed out by laplacian. to get it, assume laplasian=0 and solve phase-field eq for H
H[defect_i] = (p.gc*d[defect_i])/(2*p.l*(1.0-d[defect_i]))

F = 10.0 

u_sol, d_sol, H_sol, E_vec = main_solver(x,u,d,H,F,p)

# save all iterations
serialize(p.savefile, (x, u_sol, d_sol, H_sol, E_vec))