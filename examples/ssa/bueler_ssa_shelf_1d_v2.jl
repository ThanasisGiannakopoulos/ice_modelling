using Plots
using Printf
using Serialization

"""
Solve the SSA equation for a floating shelf, using a fixed-point (Picard-type) iteration.

- The unknown is the ice flow speed u.
- The computational domain is x=[x_grounding, x_calving] = [0,x_max]
- It is a boundary value problem with boundary data:
    - u(0) = ug
    - u_x(x_max) = gamma = [rho*g/(4B)*(1-rho/rhow)*H(x_max)]^n, with H the ice thickness.
- The algortihm needs an initial guess u0, which is u0 = gamma*x.
- For each iteration the solution is u = A^{-1}*b, where A is a matrix and b, u vectors.
"""

Base.@kwdef struct constants
    rho::Float64 = 900.0 # kg m^-3
    rhow::Float64 = 1000.0 # kg m^-3
    n::Float64 = 3.0 # for Glen's law
    B::Float64 = 1.9e8 # Pa s^1/3
    g::Float64 = 9.8 # m s^-2
    secpera::Float64 = 365*24*3600 # seconds per one year
    A::Float64 = 1.4579e-25 # same as 1/B^3
end

# midpoint, in vector N, out vector N-1
function stagav(f)
    N = length(f)
    f_mid = zeros(N-1)
    for i in 1:N-1
        f_mid[i] = 0.5*(f[i]+f[i+1])
    end
    return f_mid
end

# derivative on collocated grid
function regslope(dx,f)
    N = length(f)
    f_x = zeros(N)
    idx = 1/dx
    f_x[1] = idx*(f[2]-f[1]) # forward FD
    for i in 2:N-1
        f_x[i] = 0.5*idx*(f[i+1]-f[i-1]) # centred FD
    end
    f_x[end] = idx*(f[end]-f[end-1]) # backward FD
    return f_x
end

# derivative on staggerd grid
function stagslope(dx,f)
    N = length(f)
    f_x = zeros(N-1)
    idx = 1/dx
    for i in 1:N-1
        f_x[i] = idx*(f[i+1]-f[i]) # forward FD
    end
    return f_x
end

# build matrix A and vector b; see notes eq 29 and below.
function build_A_b(dx, H, u, ug, cst::constants)
    
    # load parameters
    n = cst.n
    rho = cst.rho
    rhow = cst.rhow
    B = cst.B
    g = cst.g
    A = cst.A

    # defs
    N = length(H)
    
    # for rhs
    r = rho/rhow
    b = -r * H
    h = H .+ b 
    beta = rho * g * H .* regslope(dx,h)
    # for BC at calving front
    gamma = (0.25 * A^(1/n) * (1 - rho/rhow) * rho * g * H[end] )^n
    Hstag = stagav(H) # length N-1
    eps_reg = (1.0 / cst.secpera) / x[end]
    uxstag = stagslope(dx, u) # length N-1
    sqr_ux_reg = uxstag.^2 .+ eps_reg^2 # length N-1
    W = zeros(N) 
    W[1:end-1] = 2 * A^(-1/n) * Hstag .* (sqr_ux_reg).^(((1/n)-1)/2)
    W[end] = copy(W[end-1])
    
    rhs = dx^2 * beta
    rhs[1] = ug
    rhs[end] = rhs[end] - 2*gamma*dx*W[end]

    # initiate A, b
    A = zeros(N,N)
    
    # set first row for BC u(0) = ug
    A[1,1] = 1.0
    
    # inner points
    for i in 2:N-1
        A[i,i-1] = W[i-1]
        A[i,i] = -(W[i-1]+W[i])
        A[i,i+1] = W[i]
    end
    A[end,end-1] = W[end-1] + W[end]
    A[end,end] = -(W[end-1] + W[end])
    
    return A, rhs
end

# function find_u(A, b)
#     return A \ b
# end
function find_u(A, b)
    # row scaling for numerical stability
    scale = vec(maximum(abs, A, dims=2))
    N = size(A,1)
    for i in 1:N
        A[i, :] .= A[i, :] ./ scale[i]
    end
    b .= b ./ scale

    # solve
    return A \ b
end

function solve_u(x, H, u0, ug, tol, max_iter, cst::constants, print_every, save_every)

    # grid spacing
    dx = x[2] - x[1]
    # make [] to hold solution and iterations
    u_list = []
    # initiate u
    u = u0
    append!(u_list, [u])
    err = 1.0e3

    # start time evolution
    iter = 0
    while err >= tol && iter < max_iter

        # calculate A, b
        A, b = build_A_b(dx, H, u, ug, cst::constants)
        # find u and store it in uf
        unew = find_u(A, b)
        # err
        err = maximum(abs.(unew .-u))
        u = unew
        iter+=1
        # print message
        if iter % print_every == 0
            @printf("it = %d , err = %.3e\n", iter, err)
        end

        # save u
        if iter % save_every == 0
            append!(u_list, [copy(u)])
        end

    end # end while
        
    return u_list, H, iter
end # end solve_u

# initial guess for u0
function init_u(x, H, cst::constants, ug)
    """
    https://github.com/bueler/karthaus/blob/master/mfiles/ssainit.m
    linear initial guess from u_xx = 0, u(0)=0, and fraction of calving condition;
    APPROPRIATE TO ICE SHELVES
    """
    n = cst.n
    rho = cst.rho
    rhow = cst.rhow
    B = cst.B
    g = cst.g
    C = g*rho*(1-rho/rhow)
    gamma = cst.A*(0.25*C*H[end])^n

    return gamma*x
end

# analytical solution
function analytical(x, cst::constants, ug, Hg, M0)
    rho = cst.rho # kg m^-3
    rhow = cst.rhow # kg m^-3
    n = cst.n # for Glen's law
    B = cst.B # Pa s^1/3
    g = cst.g # m s^-2

    # for convenience
    r = rho/rhow
    Cs = cst.A*(0.25*rho*g*(1-r))^n
    qg = ug*Hg

    # initiate u, H
    Nx = length(x)
    u = zeros(Nx)
    H = zeros(Nx)

    # form u, H
    for i in 1:Nx
        u[i] = ( ug^(n+1) + (Cs/M0) * ((M0 * x[i] + qg)^(n+1) - qg^(n+1)) )^(1/(n+1))
        H[i] = (M0 * x[i] + qg) / u[i];
    end

    return H, u
end

# grid
res = 6 # double resol
Nx = 2^res*64 
L  = 200.0e3
x = range(0, L, length=Nx)

# load constants
cst = constants()

# for analytical solution (and boundary condition at x=0)
ug = 50.0/cst.secpera # m/s; speed at grounding line
Hg = 500.0 # m; ice thickness at grounding line
M0 = 0.3/cst.secpera # m/s; surface mass balance
# analytical sollution
H, u_analytic = analytical(x, cst::constants, ug, Hg, M0)

print_every = 1#10*2^res
save_every  = 1#10*2^res

u0 = init_u(x, H, cst, ug)
tol = 1.0e-14
max_iter = 1e2
# under relaxation factor (0 < theta <= 1)
#theta = 0.5
u_list, H_used, iter = solve_u(x, H, u0, ug, tol, max_iter, cst::constants, print_every, save_every)

# save
# mkdir if it does not exist already
if isdir("./runs/")==false
    mkdir("./runs/")
end
serialize("runs/solution_$(res)res.jls", (x, u_list, H_used, iter))