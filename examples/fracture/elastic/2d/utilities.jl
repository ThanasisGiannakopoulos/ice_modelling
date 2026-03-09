# ------------------------------------------------------------
# Grid indexing (column major)
# ------------------------------------------------------------
lin(i,j,Nx) = i + (j-1)*Nx

# ------------------------------------------------------------
# Derivatives (1st order boundary, 2nd interior)
# ------------------------------------------------------------
function Dx(f,i,j,p::Params)
    Nx = p.Nx
    Δx = p.Δx
    if i == 1
        return (f[i+1,j] - f[i,j]) / Δx
    elseif i == Nx
        return (f[i,j] - f[i-1,j]) / Δx
    else
        return (f[i+1,j] - f[i-1,j]) / (2Δx)
    end
end

function Dy(f,i,j,p::Params)
    Ny = p.Ny
    Δy = p.Δy
    if j == 1
        return (f[i,j+1] - f[i,j]) / Δy
    elseif j == Ny
        return (f[i,j] - f[i,j-1]) / Δy
    else
        return (f[i,j+1] - f[i,j-1]) / (2Δy)
    end
end

function Dxx(f,i,j,p::Params)
    Nx = p.Nx
    Δx = p.Δx
    if i == 1
        return (f[1,j] - 2*f[2,j] + f[3,j]) / Δx^2
    elseif i == Nx
        return (f[Nx-2,j] - 2*f[Nx-1,j] + f[Nx,j]) / Δx^2
    else
        return (f[i-1,j] - 2*f[i,j] + f[i+1,j]) / Δx^2
    end
end

function Dyy(f,i,j,p::Params)
    Ny = p.Ny
    Δy = p.Δy
    if j == 1
        return (f[i,1] - 2*f[i,2] + f[i,3]) / Δy^2
    elseif j == Ny
        return (f[i,Ny-2] - 2*f[i,Ny-1] + f[i,Ny]) / Δy^2
    else
        return (f[i,j-1] - 2*f[i,j] + f[i,j+1]) / Δy^2
    end
end

# degradation function
function g(d,i,j,p::Params)
    return (1.0-d[i,j])^2 + p.k
end

"""
    To solve for displacement a with traction condition on left and right we need to remove rigid body motion: translation in x,y and rotation. These are 3 dofs. the matrix M for M*u=b, with u=[a1,a2] has rank=Ntot-3, with Ntot=Nx*Ny, in this setup. 

    One way to deal with this is to impose additional constraints for the mean a1, a2 displacements, as well as the mean rotation of the square.

    The idea is to numerically approximate the following:
    mean(a1), mean(a2), and mean(rotation) with rotation = x*a2 - y*a1, and set them equal to the induced translations and rotations from the boundary conditions on left/right.

    Works well for symmetric f1_left/f1_right (strech/squeeze) and gives the expected Poisson ration.

    The rest of the combinations seem off:
    e.g. if I apply the same f1 left and right (with the same direction), I don't get just translation, but I get deformation. Similary for f2. 
    TODO: Check this more.
"""
function build_constraints(p::Params)

    Nx, Ny = p.Nx, p.Ny
    N = Nx*Ny
    C = zeros(2N, 3)

    lin(i,j) = i + (j-1)*Nx

    for j in 1:Ny, i in 1:Nx
        k = lin(i,j)

        # mean ux
        C[k,1] = 1.0/N

        # mean uy
        C[N+k,2] = 1.0/N

        # rotation
        C[k,3]      = -y[j]/N
        C[N+k,3]    =  x[i]/N
    end

    return C
end