"""
Taking 2D Wave NBC, refactoring and cleaning
"""

# Useful imports
using LinearAlgebra
using Plots
using SparseArrays
# using Pkg 
using BenchmarkTools

# Local files for this project
include("./methods.jl") # Pull in Euler and RK2 Explicit Methods
include("./helper_functions.jl") # Some nice stack, unpack, etc methods
include("./sbp_ops.jl") # We'll actually build our ops
include("./plotting_2D.jl")
include("./convergence_testing.jl")

# Globals

# Wave speed
C = pi / 5

# Exact term for Manufactured Solution
function exact(t, y, z)
    # Dear god plz be merciful
    return sin(C*(y+z) - t)
end

# Exact term for Manufactured Solution Partial Derv with respect to either y or z
function exact_derv(t, y, z)
    # Dear god plz be merciful
    return C*cos(C*(y+z) - t)
end

function source_term(t, y_mesh, z_mesh)
    ny = length(y_mesh)
    nz = length(z_mesh)
    res = zeros(length(y_mesh) * length(z_mesh))
    for i in eachindex(y_mesh)
        for z in eachindex(z_mesh)
            inside = C*(y_mesh[i] + z_mesh[z]) - t
            res[(i-1)*nz + z] = (-1 * sin(inside)) - (-2*(C^4) * sin(inside))
        end
    end
    return res
end

function initialize(y_mesh, z_mesh)
    ny = length(y_mesh)
    nz = length(z_mesh)
    N = ny*nz
    res = zeros(2*N)

    # now try u(t, y, z) = t + y
    # u(t=0) = 0 + y
    # v(t=0) = 1
    for i in eachindex(y_mesh)
        for j in eachindex(z_mesh)
            inside = C*(y_mesh[i] + z_mesh[j])
            res[(i-1)*nz + j] = sin(inside) # U init
            res[(i-1)*nz + j + N] = -cos(inside )
        end
    end
    
    return  res # ST F = U_tt - c^2Uxx
end

# Boundary Conditions:

function g_y(y_mesh, z, t)
    # Displacement U(0, Z, t)
    return [exact_derv(t, y, z) for y in y_mesh]
end

function g_z(z_mesh, y, t)
    return [exact(t, y, z) for z in z_mesh]
end


## Let's Start with all Dirichlet for now
function g_prime(offset, mesh, t)
     # Velocity V(0, Z, t)
     res = zeros(length(mesh))
     for i in eachindex(mesh)
         res[i] = 2*t
     end
     return res
 end




# Setup SAT Penalty Terms
# Dirichlet is working for y=0, y=Ny so apply for z=0, z=Nz as well
function p_f(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2y, D2z) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(HIy, Iz) * Ef * (u_res[1, :] .- g_z(z_mesh, y_mesh[1],  t))
    t2 = beta .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Ef * (u_res[1, :] .- g_z(z_mesh, y_mesh[1], t))
    # return zeros(length(t1))
    return t1 .+ t2  
end

function p_r(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2y, D2z) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(HIy, Iz) * Er * (u_res[end, :] .- g_z(z_mesh, y_mesh[end],  t))
    t2 = beta .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Er * (u_res[end, :] .- g_z(z_mesh, y_mesh[end], t))
    # return zeros(length(t1))
    return t1 .+ t2  

end

function p_d(params, x, t)

    # Grab parameters from operators to vectors
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2y, D2z) = params

    # Set Constants per Erickson and Dunham 2014
    alpha_d = -1


    # Massage UV vector so this is a bit easier to work with
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)

    tmp = zeros((Ny+1) * (Nz+1) * 2) # Pack it this way to fit into existing framework


    tmp[1: (Ny+1) * (Nz+1)] = kron(Iy, BSz) * x[1: (Ny+1) * (Nz+1)]
    convert!(tmp, y_mesh, z_mesh, u_res, v_res) # pass the total product into the u_res

    # First term in the sum
    t1 = alpha_d .* kron(Iy, HIz) * Ed * (u_res[:, end] .- g_y(y_mesh, z_mesh[end], t))
    return t1  

end

function p_s(params, x, t)

     # Grab parameters from operators to vectors
     (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2y, D2z) = params

     # Set Constants per Erickson and Dunham 2014
     alpha_s = -1

 
     # Massage UV vector so this is a bit easier to work with
     u_res = zeros(Ny+1, Nz+1)
     v_res = zeros(Ny+1, Nz+1)

     tmp = zeros((Ny+1) * (Nz+1) * 2) # Pack it this way to fit into existing framework


     tmp[1: (Ny+1) * (Nz+1)] = kron(Iy, BSz) * x[1: (Ny+1) * (Nz+1)]
     convert!(tmp, y_mesh, z_mesh, u_res, v_res) # pass the total product into the u_res
 
     # First term in the sum
     t1 = alpha_s .* kron(Iy, HIz) * Es * (u_res[:, 1] .+ g_y(y_mesh, z_mesh[1], t))
    
     return t1  

end


function rhs(t, x, params)
    # Total right hand side of our ODEs
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2y, D2z) = params
    N = (Ny + 1) * (Nz + 1)
    u = x[1:N]
    v = x[N+1:2*N]
    b = zeros(2 * N)
    b[1:N] = v # move u = v part
    b[N+1:2*N] = C^2 * ((D2y*u + D2z*u)) .+ p_f(params, x, t) .+ p_r(params, x, t) + p_s(params, x, t) .+ p_d(params, x, t) .+ source_term(t, y_mesh, z_mesh)# update v with sbp
    return b
end

function run(dy, dz, dt)
    
    # MESHING

    # Space
    Y0 = 0
    Z0 = 0

    YN = 5
    ZN = 5

    DY = dy
    DZ = dz

    Y_GRID = Y0:DY:YN
    Z_GRID = Z0:DZ:ZN # No wayyyyyyy

    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    N = (NY + 1) * (NZ + 1)

    # Time
    A = 0
    B = 5
    DT = dt

    T_GRID = A:DT:B
    NT = length(T_GRID)

    # Get SBP Operators
    print("Timing for SBP OP Creation:\n")
    @time (D2y, D2z, Iy, Iz, Hy, Hz, HIy, HIz, BSy, BSz, mu, Ef, Er, Es, Ed) = sbp_operators(Y0, YN, Z0, ZN, NY, NZ, DY, DZ)
    #     (sparse(kron(D2y_test, Iz)), sparse(kron(Iy, D2z_test)), Iy, Iz, Hy_test, Hz_test, HIy_test, HIz_test, sparse(S0_y+SN_y), sparse(S0_z+SN_z), mu, Ef, Er, Es, Ed)
    # Set Up Initials
    print(size(D2y), size(D2z))
    result = zeros(2 * (NY+ 1) * (NZ + 1), NT)
    c = initialize(Y_GRID, Z_GRID)
    params = (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, DY, DZ, NY, NZ, Y_GRID, Z_GRID, D2y, D2z)

    #        (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh,D2y, D2z) = params

    #forward_euler!(f, c, dt, result, time_mesh, params)
    print("\nTiming for RK2:\n")
    @time rk2!(rhs, c, DT, result, T_GRID, params)

    plot_2d!(Y_GRID, Z_GRID, T_GRID, result, exact)

    # Only return the final time result for error
    return result[:, end]
end

converge_2D(exact)

        



