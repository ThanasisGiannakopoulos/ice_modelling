using Plots
using Printf
using Serialization


"""
Solve the Shallow Ice Approximation (SIA) with a flat Bed
∂_t h = ∂_x*(D ∂_x*h)
with D = Γ*h^5*|∂_x*h|^2
where h is evaluated at midpoints (i+1/2) and ∂_x is 1st order accurate FD
Γ is constant (includes the rate A which is assumed temperature independent here)

Analytical solution: Halfar
see Karthaus notes chap. 8
"""

# build D = Gamma h^5 * h_x^2 for d/dx*(D*d/dx*(h)) and adaptive step
function build_D(Gamma, h, dx)
    Nx = length(h)
    Drt = zeros(Nx)
    Dlt = zeros(Nx)
    one_over_dx2 = 1/dx^2
    for i in 2:Nx-1
        Hrt = 0.5*(h[i+1] + h[i]) # right staggered point from i
        Hlt = 0.5*(h[i] + h[i-1]) # left staggered point from i
        a2rt = one_over_dx2*(h[i+1] - h[i])^2 # del^2 @ right staggered from i
        a2lt = one_over_dx2*(h[i] - h[i-1])^2 # del^2 @ left staggered from i
        Drt[i] = Gamma * Hrt^5 * a2rt
        Dlt[i] = Gamma * Hlt^5 * a2lt
    end

    return Drt, Dlt
end

function rhs(Gamma, h, dx)
    Nx = length(h)
    rhs = zeros(Nx)
    one_over_dx2 = 1/dx^2
    for i in 2:Nx-1
        # rhs[i] = Gamma * (0.5 * (h[i+1]+h[i]))^5 * (one_over_dx2 * (h[i+1]-h[i])^2) * (one_over_dx2*(h[i+1]-h[i])) - Gamma * (0.5 * (h[i]+h[i-1]))^5 * (one_over_dx2 * (h[i] - h[i-1])^2) * (one_over_dx2 * (h[i] - h[i-1]))
        Hrt = 0.5*(h[i+1] + h[i]) # right staggered point from i
        Hlt = 0.5*(h[i] + h[i-1]) # left staggered point from i
        a2rt = one_over_dx2*(h[i+1] - h[i])^2 # del^2 @ right staggered from i
        a2lt = one_over_dx2*(h[i] - h[i-1])^2 # del^2 @ left staggered from i
        Drt = Gamma * Hrt^5 * a2rt
        Dlt = Gamma * Hlt^5 * a2lt
        rhs[i] = one_over_dx2*(Drt*(h[i+1] - h[i]) - Dlt*(h[i] - h[i-1]))
    end

    return rhs
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
    dt = 0.0
    f = 0.0*similar(h)
    while t<tf
        
        # calculate D, f, dt
        Drt, Dlt = build_D(Gamma, h, dx)
        D = vcat(Drt, Dlt)
        f = rhs(Gamma, h, dx)
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
            append!(t_list, t)
            append!(dt_list, dt)
            append!(f_list, [copy(f)])
        end
        
        # update h, t, step, f
        h .= h .+ dt.*f
        t = t + dt
        step += 1
        
    end # end while

    append!(h_list, [copy(h)])
    append!(t_list, t)
    append!(dt_list, dt)
    append!(f_list, [copy(f)])

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
    t = t / t0
    N = length(x)
    inside = zeros(N)
    for i in 1:N
        test = 1.0 - (abs(x[i]) / t^beta)^((n+1) / n)
        if test > 0
            inside[i] = test
        end # if
    end # for
    H = H0 .* inside.^(n / (2*n+1)) / t^alpha
    return H
end

# for the run
res = 1 # double resol
Nx = 2^res*64
L  = 1.7*1000e3
#L  = 100e3
x = range(-L, L, length=Nx)

secpera = 365*24*3600 # seconds per one year
t0 = 2*100*secpera
tf = 2*1e4*secpera

h0 = Halfar(t0, x)
# ---------------- Initial condition ----------------
#H = 200 .* exp.(-(x .- L/2).^2 ./ (15e3)^2)
#h0 = copy(H)

print_every = 200*2^res
save_every  = 400*2^res

h_list, t_list, dt_list, f_list = evol(x, h0, t0, tf, print_every, save_every)

# Save
# mkdir if it does not exist already
if isdir("./runs/")==false
    mkdir("./runs/")
end
serialize("runs/solution_$(res)res.jls", (x, h_list, t_list, dt_list, f_list))