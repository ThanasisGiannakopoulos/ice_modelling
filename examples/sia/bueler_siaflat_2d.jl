using Plots
using Printf
using Serialization


"""
Solve the Shallow Ice Approximation (SIA) with a flat Bed
∂_t h = ∂_x*(D ∂_x*h) + ∂_y*(D ∂_y*h)
with D = Γ*h^5*|∂_x*h + ∂_y*h|^2
where h is evaluated at midpoints (i+1/2) and ∂_x is 1st order accurate FD
Γ is constant (includes the rate A which is assumed temperature independent here)

Analytical solution: Halfar
see Karthaus notes chap. 8
"""

# build D = Gamma h^5 * (h_x^2 + h_y^2) for d/dx*(D*d/dx*(h)) and adaptive step
function build_D(Gamma, h, dx, dy)
    Nx = length(h[:,1]) # index i
    Ny = length(h[1,:]) # index j
    Drt = zeros(Nx, Ny) # right from i,j
    Dlt = zeros(Nx, Ny) # left from i,j
    Dup = zeros(Nx, Ny) # up from i,j
    Ddn = zeros(Nx, Ny) # down from i,j
    idx2 = 1/dx^2
    idy2 = 1/dy^2
    for i in 2:Nx-1
        for j in 2:Ny-1
            Hrt = 0.5*(h[i+1,j] + h[i,j])
            Hlt = 0.5*(h[i,j] + h[i-1,j])
            Hup = 0.5*(h[i,j+1] + h[i,j])
            Hdn = 0.5*(h[i,j] + h[i,j-1])
            # (d_x + d_y)^2 @ right i,j; uses midpoints for d_y and centered FD
            a2rt = idx2*(h[i+1,j]-h[i,j])^2 +
            0.0625*idy2*(h[i+1,j+1]+h[i,j+1]-h[i+1,j-1]-h[i,j-1])^2
            # (d_x + d_y)^2 @ left from i,j; uses midpoints for d_y and centered FD
            a2lt = idx2*(h[i,j]-h[i-1,j])^2 +
            0.0625*idy2*(h[i-1,j+1]+h[i,j+1]-h[i-1,j-1]-h[i,j-1])^2           
            # (d_x + d_y)^2 @ up i,j; uses midpoints for d_x and centered FD
            a2up = idy2*(h[i,j+1]-h[i,j])^2 +
            0.0625*idx2*(h[i+1,j+1]+h[i+1,j]-h[i-1,j+1]-h[i-1,j])^2
            # (d_x + d_y)^2 @ down from i,j; uses midpoints for d_x and centered FD
            a2dn = idy2*(h[i,j]-h[i,j-1])^2 +
            0.0625*idx2*(h[i+1,j]+h[i+1,j-1]-h[i-1,j]-h[i-1,j-1])^2
            Drt[i,j] = Gamma * Hrt^5 * a2rt
            Dlt[i,j] = Gamma * Hlt^5 * a2lt
            Dup[i,j] = Gamma * Hup^5 * a2up
            Ddn[i,j] = Gamma * Hdn^5 * a2dn
        end # for j
    end #for i

    return Drt, Dlt, Dup, Ddn
end

function rhs(Gamma, h, dx, dy)
    Nx = length(h[:,1]) # index i
    Ny = length(h[1,:]) # index j
    rhs = zeros(Nx, Ny)
    idx2 = 1/dx^2
    idy2 = 1/dy^2
    for i in 2:Nx-1
        for j in 2:Ny-1
            Hrt = 0.5*(h[i+1,j] + h[i,j])
            Hlt = 0.5*(h[i,j] + h[i-1,j])
            Hup = 0.5*(h[i,j+1] + h[i,j])
            Hdn = 0.5*(h[i,j] + h[i,j-1])
            # (d_x + d_y)^2 @ right i,j; uses midpoints for d_y and centered FD
            a2rt = idx2*(h[i+1,j]-h[i,j])^2 +
            0.0625*idy2*(h[i+1,j+1]+h[i,j+1]-h[i+1,j-1]-h[i,j-1])^2
            # (d_x + d_y)^2 @ up i,j; uses midpoints for d_x and centered FD
            a2up = idy2*(h[i,j+1]-h[i,j])^2 +
            0.0625*idx2*(h[i+1,j+1]+h[i+1,j]-h[i-1,j+1]-h[i-1,j])^2 
            # (d_x + d_y)^2 @ left from i,j; uses midpoints for d_y and centered FD
            a2lt = idx2*(h[i,j]-h[i-1,j])^2 +
            0.0625*idy2*(h[i-1,j+1]+h[i,j+1]-h[i-1,j-1]-h[i,j-1])^2
            # (d_x + d_y)^2 @ down from i,j; uses midpoints for d_x and centered FD
            a2dn = idy2*(h[i,j]-h[i,j-1])^2 +
            0.0625*idx2*(h[i+1,j]+h[i+1,j-1]-h[i-1,j]-h[i-1,j-1])^2
            Drt = Gamma * Hrt^5 * a2rt
            Dlt = Gamma * Hlt^5 * a2lt
            Dup = Gamma * Hup^5 * a2up
            Ddn = Gamma * Hdn^5 * a2dn
            rhs[i,j] = 
            idx2*(Drt*(h[i+1,j]-h[i,j])-Dlt*(h[i,j]-h[i-1,j])) +
            idy2*(Dup*(h[i,j+1]-h[i,j])-Ddn*(h[i,j]-h[i,j-1]))
        end # for j
    end #for i

    return rhs
end

function find_dt(D, dx, dy)
    maxD = maximum(D)
    dt = 0.25*min(dx,dy)^2/maxD
    return dt
end

# KO diss; 2nd order
function add_KO2!(rhs, h, dx, dy; eps = 0.01)
    Nx, Ny = size(h)
    idx2 = 1/dx
    idy2 = 1/dy
    i16 = 1/16
    for i in 3:Nx-2
        for j in 3:Ny-2
            diss =
                (h[i+2,j] - 4*h[i+1,j] + 6*h[i,j] - 4*h[i-1,j] + h[i-2,j]) * idx2 * i16
                +
                (h[i,j+2]- 4*h[i,j+1] + 6*h[i,j] - 4*h[i,j-1] + h[i,j-2]) * idy2 * i16
            rhs[i,j] += -eps * diss
        end
    end
end

# evolution; give time in seconds
function evol(x, y, h0, t0, tf, print_every, save_every)

    # constants
    secpera = 365*24*3600 # seconds per one year
    g       = 9.81
    rho     = 910.0
    A       = 1e-16/secpera # rate factor, assumed temperature independent in this SIA model
    Gamma   = 2*A*(rho*g)^3/5   # for n=3 Glenn's law

    # grid spacing
    dx = x[end] - x[end-1]
    dy = y[end] - y[end-1]
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
        Drt, Dlt, Dup, Ddn = build_D(Gamma, h, dx, dy)
        D = vcat(Drt, Dlt, Dup, Ddn)
        f = rhs(Gamma, h, dx, dy)
        #println(maximum(D))
        dt = find_dt(D, dx, dy)
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
        
        # apply 2nd order KO diss
        #add_KO2!(f, h, dx, dy; eps = 0.001)
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
function Halfar(t,x,y)
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

    t = t / t0
    Nx = length(x)
    Ny = length(y)
    inside = zeros(Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            r = sqrt(x[i]^2 + y[j]^2)
            r = r / R0
            test = 1.0 - (r / t^beta)^((n+1) / n)
            if test > 0
                inside[i,j] = test
            end # if
        end # for i
    end # for j
    H = H0 .* inside.^(n / (2*n+1)) / t^alpha
    return H
end

# for the run
res = 2 # double resol
Nx = 2^res*64 +1
Ny = 2^res*64 +1
L  = 1.2*1000e3
#L  = 100e3
x = range(-L, L, length=Nx)
y = range(-L, L, length=Ny)

secpera = 365*24*3600 # seconds per one year
t0 = 2*100*secpera
tf = 2*1e4*secpera

h0 = Halfar(t0, x, y)

print_every = 50*2^res
save_every  = 100*2^res

h_list, t_list, dt_list, f_list = evol(x, y, h0, t0, tf, print_every, save_every)

# Save
# mkdir if it does not exist already
if isdir("./runs/")==false
    mkdir("./runs/")
end
serialize("runs/2d_solution_$(res)res.jls", (x, y, h_list, t_list, dt_list, f_list))