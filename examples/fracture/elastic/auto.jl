using LinearAlgebra
using SparseArrays
using Plots

# =========================================================
# 1. Parameters
# =========================================================
L   = 1.0            # rod length
N   = 200*2            # number of FV cells
dx  = L / N

E   = 1.0            # Young's modulus
Gc  = 1e-2           # fracture toughness
ℓ   = 0.02/4           # length scale
κ   = 1e-8           # residual stiffness

nsteps = 150*2         # load steps
tol    = 1e-6
maxit  = 50

ū(t) = 2*0.5 * t     # imposed displacement at right end

# =========================================================
# 2. Grid and fields
# =========================================================
x = collect(0 : dx : L-dx)#collect(dx/2 : dx : L - dx/2)

u = zeros(N)         # displacement
d = zeros(N)         # phase field
H = zeros(N)         # history field

# -----------------------------
# Initial defect as small displacement
# -----------------------------
x0 = 0.5             # defect center
σ  = 0.05            # width of defect

for i in 1:N
    H[i] = 0.2 * exp(-0.5*((x[i]-x0)/σ)^2) # 0.1*rand()
end

# Storage for all steps
u_hist = zeros(N, nsteps)
d_hist = zeros(N, nsteps)
H_hist = zeros(N, nsteps)

ubar_hist = zeros(nsteps)
reaction_hist = zeros(nsteps)

# =========================================================
# 3. Helper functions
# =========================================================
g(d) = (1 - d)^2 + κ
εplus(ε) = max(ε, 0.0)

# =========================================================
# 4. Mechanical solver (FV)
# =========================================================
function solve_mechanics!(u, d, ubar)
    A = spzeros(N, N)
    b = zeros(N)

    for i in 1:N
        # left face
        if i == 1
            uL = -u[1]          # u(0)=0
            dL = d[1]
        else
            uL = u[i-1]
            dL = d[i-1]
        end
        εL = (u[i] - uL) / dx
        cL = E/dx^2 * (εL > 0 ? g(0.5*(d[i]+dL)) : 1.0)

        # right face
        if i == N
            uR = 2*ubar - u[N]   # u(L)=ubar
            dR = d[N]
        else
            uR = u[i+1]
            dR = d[i+1]
        end
        εR = (uR - u[i]) / dx
        cR = E/dx^2 * (εR > 0 ? g(0.5*(d[i]+dR)) : 1.0)

        # matrix assembly
        A[i,i] += cL + cR
        if i > 1
            A[i,i-1] -= cL
        end
        if i < N
            A[i,i+1] -= cR
        end

        if i == N
            b[i] += ubar#cR * 2*ubar
        end
    end
    A[1,1] = 1.0
    A[1,2] = 0.0
    A[N,N] = 1.0
    A[N,N-1] = 0.0
    #b[1] += -ubar

    u[:] = A \ b
end

# =========================================================
# 5. Phase-field solver (FV)
# =========================================================
function solve_phasefield!(d, H)
    A = spzeros(N, N)
    b = zeros(N)

    for i in 1:N
        A[i,i] = 2*Gc*ℓ/dx^2 + Gc/ℓ + 2*H[i]

        if i > 1
            A[i,i-1] = -Gc*ℓ/dx^2
        else
            A[i,i] -= Gc*ℓ/dx^2   # d'(0)=0
        end

        if i < N
            A[i,i+1] = -Gc*ℓ/dx^2
        else
            A[i,i] -= Gc*ℓ/dx^2   # d'(L)=0
        end

        b[i] = 2*H[i]
    end

    d[:] = A \ b #clamp.(A \ b, 0.0, 1.0)
end

# =========================================================
# 6. Load stepping + staggered solve
# =========================================================
for step in 1:nsteps
    t = step / nsteps
    ubar = ū(t)
    ubar_hist[step] = ubar

    for it in 1:maxit
        u_old = copy(u)
        d_old = copy(d)

        # 1. Mechanics solve
        solve_mechanics!(u, d, ubar)

        # 2. History update (positive strain energy)
        for i in 1:N
            ε = i == 1     ? (u[2]-u[1])/dx :
                i == N     ? (u[N]-u[N-1])/dx :
                             (u[i+1]-u[i-1])/(2dx)

            ψ⁺ = 0.5 * E * εplus(ε)^2
            H[i] = max(H[i], ψ⁺)
        end

        # 3. Phase-field solve
        solve_phasefield!(d, H)

        # 4. Check convergence
        err = maximum(abs.(u - u_old)) + maximum(abs.(d - d_old))
        err < tol && break
    end

    # store history
    u_hist[:,step] = u
    d_hist[:,step] = d
    H_hist[:,step] = H

    # reaction at right end
    εR = (2*ubar - 2*u[end])/dx
    reaction_hist[step] = E * (εR > 0 ? g(d[end]) * εR : εR)

    println("step $step / $nsteps   max(d) = $(maximum(d))")
end

# =========================================================
# 7. Plots
# =========================================================
# Final fields
p1 = plot(x, u, xlabel="x", ylabel="u", title="Final displacement", lw=2)
p2 = plot(x, d, xlabel="x", ylabel="d", title="Final phase field", lw=2, ylim=(0,1))
p3 = plot(x, H, xlabel="x", ylabel="H", title="Final history field", lw=2)

plot(p1, p2, p3, layout=(3,1))

# Load–displacement curve
plot(ubar_hist, reaction_hist,
     xlabel="Imposed displacement", ylabel="Reaction force",
     title="Load–displacement curve", lw=2)

# Animation of evolution
anim = @animate for step in 1:nsteps
    p1 = plot(x, u_hist[:,step], ylim=(0, maximum(u_hist)),
              ylabel="u", title="Displacement – step $step", lw=2)
    p2 = plot(x, d_hist[:,step], ylim=(0,1),
              ylabel="d", title="Phase field", lw=2)
    p3 = plot(x, H_hist[:,step], ylabel="H", title="History field", lw=2)
    plot(p1, p2, p3, layout=(3,1))
end

gif(anim, "phasefield_1d_displacement_defect.gif", fps=10)
