using LinearAlgebra, Parameters, Serialization

# Define parameters
Base.@kwdef struct Params
    left_x_BC::String = "traction" # or " Dirichlet"
    right_x_BC::String = "traction" # or " Dirichlet"
    y_BC::String = "traction_free" # "Neumann" or "traction_free"
    k::Float64 = 1e-8       # small numerical factor
    λ::Float64 = 1.0
    μ::Float64 = 1.0
    Lx::Float64 = 1.0
    Nx::Int = 15
    Δx::Float64 = Lx/(Nx-1) # spatial step
    Ly::Float64 = 1.0
    Ny::Int = 16
    Δy::Float64 = Ly/(Ny-1) # spatial step
    T::Int = 10 + 1         # number of steps for pseudotime
    dt::Float64 = 1/(T-1)   # step in pseudotime
    tol::Float64 = 1e-6     # convergence tolerance
    max_iter::Int = 20      # max staggered iterations
    savefile::String = "runs/sol.jls" # file to save iterations
    l::Float64 = 0.1 # regularization length for phase field
    C3::Float64 = 1.0 # constant for phase field eq.
end

# FD operators, index operator (matrix to list), degradation function
include("./utilities.jl")

# coefficients Q,P, that multiply ∂_x and ∂_y (respectively) of a1, a2
include("./coeffs.jl")

# returns M,b for M*a=b, where a=[a1,a2]^T is the elastic displacement in x,y, respectively
include("./displacement.jl")

# for phase field, history, etc
include("./damage.jl")


# solver for a single timestep
function ti_solver(
    a1::Matrix, a2::Matrix, 
    d::Matrix, H_old::Matrix,
    t::Float64, p::Params)

    H = copy(H_old)

    for iter in 1:p.max_iter
        # copy to do convergence check at the end
        a1_old = copy(a1)
        a2_old = copy(a2)
        d_old = copy(d)

        # get coeffs
        Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222 =
        build_QP(d_old, a1_old, a2_old,p)

        # boundary data for left/right boundaries of square domain
        f1_left, f2_left, f1_right, f2_right = displacement_BC_left_right(t, p)

        # build system for displacement
        M,b = build_displacement_system(
            Q111,P111,Q112,P112,
            Q121,P121,Q122,P122,
            Q221,P221,Q222,P222,
            p,
            f1_left,f1_right,
            f2_left,f2_right
            )

        # solve for a
        if p.left_x_BC=="traction" && p.right_x_BC=="traction"
            C = build_constraints(p)
            gg = zeros(3)
            # induced aggregate x translation form x_BC
            gg[1] = 0.0#0.5*(sum(f1_left) + sum(f1_right))/p.Ny
            # induced aggregate y translation form x_BC
            gg[2] = 0.0#0.5*(sum(f2_left) + sum(f2_right))/p.Ny
            # induced aggregate rotation from x_BC: x*a2 - y*a1 is rotation
            gg[3] = 0.0#0.5*(sum(x[1]*f2_left .- y.*f1_left) + sum(x[end]*f2_right .- y.*f1_right))/p.Ny
            K = [M  C;
                C' zeros(3,3)]
            rhs = [b;
                gg]
            a_sol = K \ rhs
        else
            a_sol = M \ b
        end
        Ntot = p.Nx*p.Ny
        # make a1,a2 sols in matrices (x,y, grid functions)
        a1 = reshape(a_sol[1:Ntot], p.Nx, p.Ny)
        a2 = reshape(a_sol[Ntot+1:2*Ntot], p.Nx, p.Ny)

        # get history
        H = history(H_old, a1, a2, p)
        println("   Iteration $iter: max(d_old) = ", maximum(d_old))

        # phase-field
        Md, bd = build_phase_field_sys(H, p)
        d_sol = Md \ bd
        d = reshape(d_sol[1:Ntot], p.Nx, p.Ny)

        # convergence check
        res = maximum(abs.(a1 .- a1_old)) + maximum(abs.(a2 .- a2_old)) + + maximum(abs.(d .- d_old))
        if res < p.tol
            println("   res = $res | converged at iteration $iter")
            break
        end # end if

    end # end iters

    return a1, a2, d, H
end # end ti_solver


# apply ti_solver at each ti in the pseudotime domain
function pseudotime_solver(
    a1_0::Matrix, a2_0::Matrix, 
    d0::Matrix, H0::Matrix, 
    p::Params)

    Nx, Ny = p.Nx, p.Ny # number of spatial nodes
    T = p.T # number of pseudotemporal steps

    # pseudotime
    t = range(p.dt,1,T-1)
    # initiale matrices to save data
    a1_all = zeros(Nx, Ny,T-1)
    a2_all = zeros(Nx, Ny,T-1)
    d_all = zeros(Nx, Ny,T-1)
    H_all = zeros(Nx, Ny,T-1)

    # save ID
    a1_all[:,:,1] = a1_0
    a2_all[:,:,1] = a2_0
    d_all[:,:,1] = d0
    H_all[:,:,1] = H0
    
    # prepare to first steps ti=1 that is t=0
    a1 = a1_0
    a2 = a2_0
    d = d0
    H_old = H0

    for ti in 2:T-1
        # get time
        tt = t[ti]
        println("ti = $ti | t = $tt:")
        # solve for ti
        a1, a2, d, H = ti_solver(a1, a2, d, H_old, tt, p)

        # save output at the end of iterations
        a1_all[:,:,ti] = copy(a1)
        a2_all[:,:,ti] = copy(a2)
        d_all[:,:,ti] = copy(d)
        H_all[:,:,ti] = copy(H)
        
        # update history for the next step
        H_old = copy(H)
    end

    return a1_all, a2_all, d_all, H_all
end


# Run example
p = Params(    
    left_x_BC = "traction", # "traction" 
    right_x_BC = "traction", # or "traction"
    y_BC="traction_free",
    λ = 1.0,
    μ = 1.0,
    Nx = 30,
    Ny = 30,
    max_iter = 10,
    T = 40 + 1,
    savefile="runs/sol_left_traction_right_traction.jls"
)

x = range(0,p.Lx,p.Nx)
y = range(0,p.Ly,p.Ny)

# ID
# displacements
a1_0 = zeros(p.Nx, p.Ny)
a2_0 = zeros(p.Nx, p.Ny)

# phase field Gaussian*step_function
d0 = zeros(p.Nx, p.Ny)
f0 = 0.5*(1.0 .+ tanh.(50*(y .- 0.5)))
g0 = 1.0*exp.(-(x .- 0.5).^2/(2*0.005))
for j in 1:p.Ny
    for i in 1:p.Nx
        d0[i,j] = g0[i]*f0[j]
    end
end

# history compatible with damage
H0 = history_from_d(d0, p)

a1_all, a2_all, d_all, H_all = pseudotime_solver(a1_0, a2_0, d0, H0, p)

# save all iterations
serialize(p.savefile, (x, y, a1_all, a2_all, d_all, H_all))
