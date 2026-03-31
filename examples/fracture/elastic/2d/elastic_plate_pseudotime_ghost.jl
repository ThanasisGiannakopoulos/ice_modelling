using LinearAlgebra, Parameters, Serialization, Dates

if isdir("./runs/")==false
    mkdir("./runs/")
end

# Define parameters
Base.@kwdef struct Params
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
    max_iter::Int = 10      # max staggered iterations
    savefile::String = "runs/sol.jls" # file to save iterations
    l::Float64 = 0.1 # regularization length for phase field
    C3::Float64 = 1.0 # constant for phase field eq.
    """ 
    parameter for star-convex energy split
    with γstar>=-1, γstar=-1 is standard model and γstar=0 -s volumetric-deviatoric split
    """
    γstar::Float64 = -1.0
    simple_d::Bool = false # solve the phase field eq w/ ∂^2 (false) or drop ∂^2 (true) 
end

# FD operators, index operator (matrix to list), degradation function
include("./utilities.jl")

# coefficients Q,P, that multiply ∂_x and ∂_y (respectively) of a1, a2
include("./coeffs.jl")

# returns M,b for M*a=b, where a=[a1,a2]^T is the elastic displacement in x,y, respectively
include("./displacement_v3.jl")

# for phase field, history, etc
include("./damage.jl")


# solver for a single timestep
function ti_solver(
    a1::Matrix, a2::Matrix, 
    d::Matrix, H_old::Matrix,
    t_past::Float64,
    t::Float64, p::Params)

    H = copy(H_old)

    for iter in 1:p.max_iter
        # copy to do convergence check at the end
        a1_old = copy(a1)
        a2_old = copy(a2)
        d_old = copy(d)

        # boundary data for left/right boundaries of square domain
        # for previous timestep
        f1_left_past, f2_left_past, f1_right_past, f2_right_past = displacement_BC_left_right(t_past, p)
        
        # # get coeffs
        # Q111,P111,Q112,P112,
        # Q121,P121,Q122,P122,
        # Q221,P221,Q222,P222 =
        # build_QP_ghost(d_old, a1_old, a2_old,p,
        # f1_left_past, f1_right_past,
        # f2_left_past, f2_right_past)
        
        # get coeffs
        Q111,P111,Q112,P112,
        Q121,P121,Q122,P122,
        Q221,P221,Q222,P222 =
        build_QP_gammastar(d_old, a1_old, a2_old,p)

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


        a_sol = M \ b
        Ntot = p.Nx*p.Ny
        # make a1,a2 sols in matrices (x,y, grid functions)
        a1 = reshape(a_sol[1:Ntot], p.Nx, p.Ny)
        a2 = reshape(a_sol[Ntot+1:2*Ntot], p.Nx, p.Ny)

        # get history
        #H = history_ghost(H_old, a1, a2, p,f1_left,f1_right)
        H = history_gammastar(H_old, a1, a2, p)
        println("   Iteration $iter: max(d_old) = ", maximum(d_old))

        # phase-field
        if p.simple_d==false
            Md, bd = build_phase_field_sys(H, p)
            d_sol = Md \ bd
            d = reshape(d_sol[1:Ntot], p.Nx, p.Ny)
        else
            d = 2.0*p.C3*p.l*H ./ (1.0 .+ 2.0*p.C3*p.l*H)
        end

        # convergence check
        res = maximum(abs.(a1 .- a1_old)) + maximum(abs.(a2 .- a2_old)) + + maximum(abs.(d .- d_old))
        println("res = $res")
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

    start_time = now()

    Nx, Ny = p.Nx, p.Ny # number of spatial nodes
    T = p.T # number of pseudotemporal steps

    # pseudotime
    t = range(0.0,1,T)#range(p.dt,1,T-1)
    # initiale matrices to save data
    a1_all = zeros(Nx, Ny,T)
    a2_all = zeros(Nx, Ny,T)
    d_all = zeros(Nx, Ny,T)
    H_all = zeros(Nx, Ny,T)

    if isfile(p.savefile)
        x_loaded, y_loaded, a1_loaded, a2_loaded, d_loaded, H_loaded =
        deserialize(p.savefile)
        last_iter = size(a1_loaded,3)
        first_iter = last_iter + 1
        # transfer loaded data to a1_all etc
        # a1_all[:,:,1:first_iter-1] = a1_loaded[:,:,1:first_iter-1]
        # a2_all[:,:,1:first_iter-1] = a2_loaded[:,:,1:first_iter-1]
        # d_all[:,:,1:first_iter-1] = d_loaded[:,:,1:first_iter-1]
        # H_all[:,:,1:first_iter-1] = H_loaded[:,:,1:first_iter-1]
        a1_all[:,:,1:last_iter] = a1_loaded
        a2_all[:,:,1:last_iter] = a2_loaded
        d_all[:,:,1:last_iter] = d_loaded
        H_all[:,:,1:last_iter] = H_loaded

        # # prepare to first step ti=first_iter-1 (last timestep the checkpoint has saved)
        # a1 = a1_loaded[:,:,end]
        # a2 = a2_loaded[:,:,end]
        # d = d_loaded[:,:,end]
    
        # H_old = H_loaded[:,:,end]
        # restart state
        a1 = a1_loaded[:,:,last_iter]
        a2 = a2_loaded[:,:,last_iter]
        d  = d_loaded[:,:,last_iter]
        H_old = H_loaded[:,:,last_iter]
        
    else

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

        first_iter = 2
    end
    
    killfile = "runs/kill"

    for ti in first_iter:T#first_iter:T-1
        
        # -------------------------
        # KILL SWITCH CHECK
        # -------------------------
        if isfile(killfile)
            println("Kill file detected. Saving solution and stopping simulation.")

            serialize(p.savefile,
                (x, y,
                 a1_all[:,:,1:ti-1],
                 a2_all[:,:,1:ti-1],
                 d_all[:,:,1:ti-1],
                 H_all[:,:,1:ti-1]))

            return a1_all[:,:,1:ti-1], a2_all[:,:,1:ti-1], d_all[:,:,1:ti-1], H_all[:,:,1:ti-1]
        end
        
        # get time
        tt_past = t[ti-1]
        tt = t[ti]
        println("ti = $ti | t = $tt:")
        # solve for ti
        a1, a2, d, H = ti_solver(a1, a2, d, H_old, tt_past, tt, p)

        # save output at the end of iterations
        a1_all[:,:,ti] = copy(a1)
        a2_all[:,:,ti] = copy(a2)
        d_all[:,:,ti] = copy(d)
        H_all[:,:,ti] = copy(H)
        
        # update history for the next step
        H_old = copy(H)
    end

    end_time = now()
    # Compute total elapsed time in minutes
    elapsed_ms = Dates.value(end_time - start_time)  # milliseconds as integer
    elapsed_minutes = elapsed_ms / 1000 / 60        # convert to minutes
    println("Total time: ", elapsed_minutes, " minutes")

    return a1_all, a2_all, d_all, H_all
end

function displacement_BC_left_right(t::Float64, p::Params)
    f1_left  = -0.0*t*1e-2*ones(p.Ny) # 0.0*ones(p.Ny)
    f2_left  = -0*1e-1*ones(p.Ny)
    
    f1_right =  0.0*t*1e-2*ones(p.Ny) # t*0.5*1e-1*ones(p.Ny)
    f2_right =  -1.0*t*1e-2*ones(p.Ny) # 0*1e-1*ones(p.Ny)
    return f1_left, f2_left, f1_right, f2_right
end

# Run example
p = Params(    
    λ = 121.15,#1.0,
    μ = 80.77,#1.0,
    Nx = 40,
    Ny = 40,
    max_iter = 10,
    T = 1*500 + 1,
    γstar = 0.0,
    k = 1e-6,
    C3 = 0.37*1e3,
    l = 0.075,
    savefile="runs/sol_gammastar_0_lambda_121.15_mu_80.77_C3_0.37*1e3_Nx40_Ny40_l0.075_ID_notch_f2right_m0.01_steps_500_simple_d_false.jls",
    simple_d = false
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
g0 = 0.95*exp.(-(x .- 0.5).^2/(2*0.00125))
# for j in 1:p.Ny
#     for i in 1:p.Nx
#         d0[i,j] = g0[i]*f0[j]
#     end
# end

# notch ID
for j in 1:p.Ny
    for i in 1:p.Nx
        #d0[i,j] = g0[i]*f0[j]
        if i == Int(p.Nx/2)
            d0[i,j] = f0[j]
        end 
    end
end

# history compatible with damage
f1_left_t0, f2_left_t0, f1_right_t0, f2_right_t0 = displacement_BC_left_right(0.0, p)
H0 = history_from_d(d0, p, f1_left_t0, f1_right_t0)

killfile = "runs/kill"

# Remove old kill file if it exists
if isfile(killfile)
    rm(killfile)
end

a1_all, a2_all, d_all, H_all = pseudotime_solver(a1_0, a2_0, d0, H0, p)

# save all iterations
serialize(p.savefile, (x, y, a1_all, a2_all, d_all, H_all))
