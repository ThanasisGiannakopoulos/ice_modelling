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
    B::Float64 = 1.9000167139607197e8#1.9e8 # Pa s^1/3
    g::Float64 = 9.8 # m s^-2
    secpera::Float64 = 365*24*3600 # seconds per one year
end

# build matrix A and vector b; see notes eq 29 and below.
function build_A_b(dx, H, u, ug, cst::constants)
    
    # load parameters
    n = cst.n
    rho = cst.rho
    rhow = cst.rhow
    B = cst.B
    g = cst.g
    
    # defs
    Nx = length(H)
    idx = 1/dx
    pow = 1/n-1 # power used in calculation of W
    C = g*rho*(1-rho/rhow)

    # initiate the matrix
    A = zeros(Nx,Nx)
    b = zeros(Nx)

    # set first row for BC u(0) = ug
    A[1,1] = 1
    b[1] = ug

    # inner points
    for i in 2:Nx-1
        # build W_{i-1/2}; left of i
        W_lt = B*(H[i]+H[i-1])*abs(idx*(u[i]-u[i-1]))^pow
        # build W_{i+1/2}; right of i
        W_rt = B*(H[i+1]+H[i])*abs(idx*(u[i+1]-u[i]))^pow
        # replace A[i,i-1], A[i,i], and A[i,i+1]
        A[i,i-1] = W_lt
        A[i,i] = -W_lt-W_rt
        A[i,i+1] = W_rt
        # build beta_i; rhs of eq
        # centred FD
        H_x = 0.5*idx*(H[i+1]-H[i-1])
        beta = 0.5*C*H[i]*H_x
        # replace b[i]
        b[i] = beta*dx^2
    end

    # last point i=N; assume W_rt = W_lt, for i=N+1/2
    # build W_{N-1/2}; left of =N
    W_lt = B*(H[end]+H[end-1])*abs(idx*(u[end]-u[end-1]))^pow
    # matrix A
    A[end,end-1] = 2.0*W_lt
    A[end,end] = -2.0*W_lt
    # vector b
    gamma = (0.25*C*H[end]/B)^n
    # backward FD at i=N
    H_x = idx*(H[end]-H[end-1])
    beta = 0.5*C*H[end]*H_x
    b[end] = -2.0*dx*W_lt*gamma + beta*dx^2

    return A, b
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

function solve_u(x, H, u0, ug, error_tol, max_iter, cst::constants, print_every, save_every)

    # constants
    secpera = cst.secpera # seconds per one year
    # make ug from m/y to m/s
    ug = ug/secpera
    # grid spacing
    dx = x[2] - x[1]
    # make [] to hold solution and iterations
    u_list = []
    # initiate ui, uf for each iteration
    ui = u0
    uf = 0.0*u0
    append!(u_list, [ui])
    # calculate error
    err = maximum(abs.(ui.-uf))

    # start time evolution
    iter = 0
    while err >= error_tol && iter < max_iter

        # calculate A, b
        A, b = build_A_b(dx, H, ui, ug, cst::constants)
        # find u and store it in uf
        uf = find_u(A, b)
        # # under-relaxation  
        # uf = theta*uf + (1-theta)*ui
        # # calculate error
        err = maximum(abs.(ui.-uf))
        # update iter
        iter+=1
        # prepare for next step by copying uf to ui
        ui = copy(uf)

        # print message
        if iter % print_every == 0
            #@printf("it = %d , err = %.14f\n", iter, err)
            @printf("it = %d , err = %.3e\n", iter, err)
        end

        # save u
        if iter % save_every == 0
            append!(u_list, [copy(uf)])
        end

    end # end while
        
    # # print last iter
    # @printf("it = %d , err = %.3e\n", iter, err)
    # # save u last iter
    # if iter % save_every == 0
    #     append!(u_list, [copy(uf)])
    # end

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
    gamma = (0.25*C*H[end]/B)^n

    return gamma*x #.+ ug/cst.secpera
end

# analytical solution
function analytical(x, cst::constants, ug, Hg, M0)
    rho = cst.rho # kg m^-3
    rhow = cst.rhow # kg m^-3
    n = cst.n # for Glen's law
    B = cst.B # Pa s^1/3
    g = cst.g # m s^-2
    secpera = cst.secpera
    ug = ug/secpera # m/s; speed at grouding line
    M0 = M0/secpera # m/s; surface mass balance

    # for convenience
    r = rho/rhow
    Cs = (0.25*rho*g*(1-r)/B)^n
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
res = 2 # double resol
Nx = 2^res*1024 
L  = 2.0e5
x = range(0, L, length=Nx)

# load constants
cst = constants()

# for analytical solution (and boundary condition at x=0)
ug = 50 # m/y; speed at grounding line
Hg = 500 # m; ice thickness at grounding line
M0 = 0.3 # m/y; surface mass balance
# analytical sollution
H, u_analytic = analytical(x, cst::constants, ug, Hg, M0)

print_every = 1#10*2^res
save_every  = 1#10*2^res

u0 = init_u(x, H, cst, ug)
error_tol = 1.0e-14
max_iter = 1e2
# under relaxation factor (0 < theta <= 1)
#theta = 0.5
u_list, H_used, iter = solve_u(x, H, u0, ug, error_tol, max_iter, cst::constants, print_every, save_every)

# save
# mkdir if it does not exist already
if isdir("./runs/")==false
    mkdir("./runs/")
end
serialize("runs/solution_$(res)res.jls", (x, u_list, H_used, iter))