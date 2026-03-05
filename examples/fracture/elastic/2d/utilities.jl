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