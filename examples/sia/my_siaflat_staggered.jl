using Plots
using Printf
using Serialization


"""
Solve the Shallow Ice Approximation (SIA) with a flat Bed
∂_t h = ∂_x*(D ∂_x*h)
with D = Γ*h^5*|∂_x*h|^2
where h is evaluated at midpoints (i+1/2) and ∂_x is 1st order accurate FD
Γ is constant (includes the rate A which is assumed temperature independent here)

Analyitical solution: Halfar
see Karthaus notes chap. 8
"""

# calculate values at i+1/2 using a linear relation
function my_midpoint(f)
    N = length(f)
    f_mid = zeros(N+1)
    #f_mid[1] = f[1]
    for i in 2:N
        f_mid[i] = 0.5 * (f[i-1] + f[i])
    end
    #f_mid[end] = f[end]
    return f_mid
end

# 1st order accurate forward finite difference
function FD1_forward(f, dx)
    N = length(f)
    f_x = zeros(N + 1)
    for i in 2:N
        f_x[i] = (f[i] - f[i-1]) / dx
    end

    #f_x[end] = 0.0
    return f_x
end

# build D = Gamma h^5 * h_x^2 for d/dx*(D*d/dx*(h)) and adaptive step
function build_D(Gamma, h, dx)
    D = zeros(length(h)+1)
    h_midpoints = my_midpoint(h)
    h_x = FD1_forward(h,dx)
    D .= Gamma * h_midpoints.^5 .* h_x.^2
    return D
end

function compute_flux(Gamma, h, dx)
    D = build_D(Gamma, h, dx)
    q = D .* FD1_forward(h,dx)
    return D, q
end

# build the rhs
# function rhs(Gamma, h, dx)
#     f = zeros(length(h))
#     D = build_D(Gamma, h, dx)
#     h_x = FD1_forward(h,dx)
#     f = FD1_forward(D.*h_x, dx)
#     # set the rhs to 0 at i=0,end (no ice flux outside domain)
#     f[1] = 0.0
#     f[end] = 0.0
#     return D, f # export D as well to use for dt adaptive
# end
function rhs(q, dx)
    Nx = length(q)-1
    f = zeros(Nx)
    for i in 1:Nx
        f[i] = (q[i+1] - q[i])/dx
    end

    return f 
end

function find_dt(D, dx)
    maxD = maximum(D)
    dt = 0.25*dx^2/maxD
    return dt
end

# evolution; give time in seconds
function evol(x, h0, t0, tf, print_every, save_every)

    # constants
    secpera = 365*24*3600 # seconds per one year
    g       = 9.81
    rho     = 910.0
    A       = 1e-16/secpera # rate factor, assumed temperature independent in this SIA model
    Gamma   = 2*A*(rho*g)^3/5   # for n=3 Glenn's law

    # grid spacing
    dx = x[2] - x[1]
    # copy h0, t0
    h = copy(h0)
    t = copy(t0)
    # make [] to hold solution and timesteps, etc
    h_list = []
    t_list = []
    dt_list = [] # adaptive timesteps
    f_list = [] # rhs

     # start time evolution
    step = 0

    while t<tf
        
        # calculate D, f, dt
        D, q = compute_flux(Gamma, h, dx)
        f = rhs(q,dx)
        #println(maximum(D))
        dt = find_dt(D, dx)
        # to not go beyond tf
        if t + dt > tf
            dt = tf - t
        end

        # print message
        if step % print_every == 0
            @printf("t = %.2f kyr, dt = %.6f yr, max(h)=%.1f m\n", t/(1000*secpera), dt/secpera, maximum(h))
        end

        # save h, t, etc
        if step % save_every == 0
            append!(h_list, [copy(h)])
            append!(t_list, t/secpera)
            append!(dt_list, dt/secpera)
            append!(f_list, [copy(f)])
        end
        
        # update h, t, step, f
        h .= h .+ dt.*f
        t = t + dt
        step+=1
        
    end # end while

    return h_list, t_list, dt_list, f_list

end # end evol

# analytical solution
function Halfar(t,x)
    # copied from https://github.com/bueler/karthaus/blob/master/mfiles/halfar.m
    g = 9.81     # constants in SI units
    rho = 910.0
    secpera = 365*24*3600 # seconds per one year
    n = 3
    A = 1.0e-16/secpera
    Gamma  = 2.0 * A * (rho * g)^3 / 5.0 # see Bueler et al (2005)

    H0 = 3600.0
    R0 = 750.0e3
    alpha = 1/9#0.1111111111111111
    beta = 1/18#0.05555555555555555
    # for t0, see equation (9) in Bueler et al (2005); result is 422.45 a:
    t0 = (beta/Gamma) * (7/4)^3 * (R0^4/H0^7)

    #r = sqrt(x.*x + y.*y)
    #r = r / R0
    # retrict in x only
    x = x / R0
    t = t*secpera / t0
    inside = zeros(length(x))
    for i in 1:length(x)
        test = 1.0 - (abs(x[i]) / t^beta)^((n+1) / n)
        if test > 0
            inside[i] = test
        end # if
    end # for
    H = H0 .* inside.^(n / (2*n+1)) / t^alpha
    return H
end

# for the run
Nx = 300
L  = 100e3
x = range(0, L, length=Nx)

secpera = 365*24*3600 # seconds per one year
t0 = 0*secpera
tf = 100000*secpera

#h0 = Halfar(t0, x)
# ---------------- Initial condition ----------------
H = 200 .* exp.(-(x .- L/2).^2 ./ (15e3)^2)
h0 = copy(H)

print_every = 100
save_every  = 200

h_list, t_list, dt_list, f_list = evol(x, h0, t0, tf, print_every, save_every)

# Save
serialize("solution.jls", (x, h_list, t_list, dt_list, f_list))