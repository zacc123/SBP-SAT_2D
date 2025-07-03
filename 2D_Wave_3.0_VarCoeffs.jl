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
include("./plotting_2D.jl") # Moved all plots to here
include("./convergence_testing.jl") # Moved Convergence here

include("./get_ops_draft_1.1.jl") # Add in the get ops now
include("./coordinate_transform.jl") # Add in coordinate tranforms

# Globals

# Wave speed
C = pi / 5

# Exact term for Manufactured Solution
function exact(t, y, z)
    # Dear god plz be merciful
    inside = yf_r(y,z) + zf_s(y,z)

    return sin(C*inside - t)
end

# Exact term for Manufactured Solution Partial Derv with respect to either y or z
function exact_derv(t, y, z)
    # Dear god plz be merciful
    inside = yf_r(y,z) + zf_s(y,z)
    return C*cos(C*inside - t)
end

function source_term(t, y_mesh, z_mesh)
    
    ny = length(y_mesh)
    nz = length(z_mesh)
    res = zeros(length(y_mesh) * length(z_mesh))
    for i in eachindex(y_mesh)
        for z in eachindex(z_mesh)
            inside = (C * (yf_r(y_mesh[i],z_mesh[z]) + zf_s(y_mesh[i],z_mesh[z]))) - t
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
            inside = C * (yf_r(y_mesh[i],z_mesh[j]) + zf_s(y_mesh[i],z_mesh[j]))
            res[(i-1)*nz + j] = sin(inside) # U init
            res[(i-1)*nz + j + N] = -cos(inside)
        end
    end
    
    return  res # ST F = U_tt - c^2Uxx
end

# Boundary Conditions:

function g_y(y_mesh, z, t)
    # Displacement U(0, Z, t)
    return [exact(t, y, z) for y in y_mesh]
end

function g_y_n(y_mesh, z, t)
    # Displacement U(0, Z, t)
    return [exact_derv(t, y, z) for y in y_mesh]
end

function g_z(z_mesh, y, t)
    return [exact(t, y, z) for z in z_mesh]
end



# Setup SAT Penalty Terms
# Dirichlet is working for y=0, y=Ny so apply for z=0, z=Nz as well
function p_f(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(HIy, Iz) * Ef * (sJ[1] .* (u_res[1, :] .- g_z(z_mesh, y_mesh[1],  t)))
    t2 = (beta) .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Ef * (sJ[1] .* (u_res[1, :] .- g_z(z_mesh, y_mesh[1], t) ))
    # return zeros(length(t1))
    return (t1 .+ t2)  
end

function p_r(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(HIy, Iz) * Er * (sJ[1] .* (u_res[end, :] .- g_z(z_mesh, y_mesh[end],  t)))
    t2 = beta .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Er *  (sJ[1] .* (u_res[end, :] .- g_z(z_mesh, y_mesh[end], t)))
    # return zeros(length(t1))
    return t1 .+ t2  

end


function p_d(params, x, t)

    # Grab parameters from operators to vectors
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params

    # Set Constants per Erickson and Dunham 2014
    alpha_d = -1


    # Massage UV vector so this is a bit easier to work with
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)

    tmp = zeros((Ny+1) * (Nz+1) * 2) # Pack it this way to fit into existing framework


    tmp[1: (Ny+1) * (Nz+1)] =  kron(Iy, BSz) * x[1: (Ny+1) * (Nz+1)]
    convert!(tmp, y_mesh, z_mesh, u_res, v_res) # pass the total product into the u_res

    # First term in the sum
    t1 = alpha_d .* kron(Iy, HIz) * Ed * (2.0 .* (0.5 .* u_res[:, end] .-  g_y_n(y_mesh, z_mesh[end], t))) 
    return  t1

end

#=
function p_d(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(Iy, HIz) * Ed * (sJ[1] .* (u_res[:, end] .- g_y(y_mesh, z_mesh[end],  t)))
    t2 = beta .* kron(Iy, HIz) * transpose(kron(Iy, BSz)) * Ed *  (sJ[1] .* (u_res[:, end] .- g_y(y_mesh, z_mesh[end], t)))
    # return zeros(length(t1))
    return t1 .+ t2  

end
=#
#=
function p_s(params, x, t)

     # Grab parameters from operators to vectors
     (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params

     # Set Constants per Erickson and Dunham 2014
     alpha_s = -1

 
     # Massage UV vector so this is a bit easier to work with
     u_res = zeros(Ny+1, Nz+1)
     v_res = zeros(Ny+1, Nz+1)

     tmp = zeros((Ny+1) * (Nz+1) * 2) # Pack it this way to fit into existing framework


     tmp[1: (Ny+1) * (Nz+1)] = (kron(Iy, BSz) * x[1: (Ny+1) * (Nz+1)])
     convert!(tmp, y_mesh, z_mesh, u_res, v_res) # pass the total product into the u_res
 
     # First term in the sum
     t1 = alpha_s .* kron(Iy, HIz) * Es * ((sJ[1] .* u_res[:, 1]) .+ g_y(y_mesh, z_mesh[1], t))
    
     return t1  

end
=#

function p_s(params, x, t)
    # Dirichlet SAT for y=Ny
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    t1 = alpha_r .* kron(Iy, HIz) * Es * (sJ[3] .* (u_res[:, 1] .- g_y(y_mesh, z_mesh[1],  t)))
    t2 = beta .* kron(Iy, HIz) * transpose(kron(Iy, BSz)) * Es *  (sJ[3] .* (u_res[:, 1] .- g_y(y_mesh, z_mesh[1], t)))
    # return zeros(length(t1))
    return t1 .+ t2  

end


function rhs(t, x, params)
    # Total right hand side of our ODEs
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    N = (Ny + 1) * (Nz + 1)
    u = x[1:N]
    v = x[N+1:2*N]
    b = zeros(2 * N)
    b[1:N] = v # move u = v part
    D = D2rr + D2sr + D2rs + D2ss
    b[N+1:2*N] = C^2 * ((D ./ J) * u) .+ p_f(params, x, t) .+ p_r(params, x, t) + p_s(params, x, t) .+ p_d(params, x, t) .+ source_term(t, y_mesh, z_mesh)# update v with sbp
    return b
end


"""
    SAT_Terms(params, x, t)

SAT Terms for SBP-SAT Formulation of 2D Wave Equation from Erickson and Dunham, et al. (2014)

Currently has 3 Dirichlet BCs, a,f,d and one Neuman s.

INPUTS:
    - params: (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ)
    - x: stacked u, v vector at current time step
    - t: current time

OUTPUT:
    - Combine stacked vector [0, ..., v] for next time step of all SAT terms.
"""
function SAT_Terms(params, x, t)
    (mu, HIy, HIz, Iy, Iz, Ef, Er, Es, Ed, BSy, BSz, dy, dz, Ny, Nz, y_mesh, z_mesh, D2rr, D2ss, D2rs, D2sr, J, sJ) = params
    
    # Dirichlet Terms
    alpha_r = -13 / dy # params from E + D 2014
    beta = 1
    
    u_res = zeros(Ny+1, Nz+1)
    v_res = zeros(Ny+1, Nz+1)
    convert!(x, y_mesh, z_mesh, u_res, v_res)

    # fault (y=y1 case)
    t1_f = alpha_r .* kron(HIy, Iz) * Ef * (sJ[1] .* (u_res[1, :] .- g_z(z_mesh, y_mesh[1],  t)))
    sat_f = t1_f .+ ((beta) .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Ef * (sJ[1] .* (u_res[1, :] .- g_z(z_mesh, y_mesh[1], t) )))

    # remote (y=Yn case)
    t1_r = alpha_r .* kron(HIy, Iz) * Er * (sJ[1] .* (u_res[end, :] .- g_z(z_mesh, y_mesh[end],  t)))
    sat_r = t1_f .+ (beta .* kron(HIy, Iz) * transpose(kron(BSy, Iz)) * Er *  (sJ[1] .* (u_res[end, :] .- g_z(z_mesh, y_mesh[end], t))))
    


end

function run_logical(rc, sc, tc, metrics, D)
    # Run the SBP SAT with Euler Post Coordinate Tranform
    # INPUTS:
        # Grid spacing info
        # rc = (r0, rN, dr)
        # sc = (s0, sN, ds)
        # tc = (t0, tN, dt)
        # metrics from the coordinate transform metrics


    # Unpack and get meshes in order
    R0, RN, DR = rc
    S0, SN, DS = sc
    T0, TN, DT = tc

    R_GRID = R0:DR:RN
    S_GRID = S0:DS:SN

    NR = length(R_GRID) - 1
    NS = length(S_GRID) - 1
    N = (NR + 1) * (NS + 1)

    T_GRID = T0:DT:TN
    NT = length(T_GRID)

    # Get SBP Operators assuming orthonormal set up on R and S
    print("Timing for SBP OP Creation:\n")
    @time (D2r, D2s, Ir, Is, Hr, Hs, HIr, HIs, BSr, BSs, mu, Ef, Er, Es, Ed) = sbp_operators(R0, RN, S0, SN, NR, NS, DR, DS)
    
    # Get Jacobians and all into correct format
    J = zeros(N)
    stack!(J, metrics.J) # Stack J into same format as R, S
    
    #->
    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    result = zeros(2 * (NR+ 1) * (NS + 1), NT)
    c = initialize(R_GRID, S_GRID)

    params = (mu, HIr, HIs, Ir, Is, Ef, Er, Es, Ed, BSr, BSs, DR, DS, NR, NS, R_GRID, S_GRID, D2rr, D2ss, D2rs, D2sr, J, metrics.sJ)


    #forward_euler!(f, c, dt, result, time_mesh, params)
    print("\nTiming for RK2:\n")
    @time rk2!(rhs, c, DT, result, T_GRID, params)

    plot_2d!(R_GRID, S_GRID, T_GRID, result, exact)

    # Only return the final time result for error
    return result[:, end]
end

function main()
    #=
    Main is the entire simulation run, starting with physical definition then to logical

    Y and Z are physical, R and S are logical coordinates
    =#
    # Define Physical Meshes
    Y0, YN, dy = (-2, 2, 0.2)
    Z0, ZN, dz = (-1, 1, 0.1)
    T0, TN, dt = (0, 2, 1e-4)

    Y_GRID = Y0:dy:YN
    Z_GRID = Z0:dz:ZN
    T_GRID = T0:dt:TN

    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    NT = length(T_GRID)

    # Now make the coordinate transform
    
    metrics = create_metrics_BP6(4, NY, NZ, yf_1, zf_1) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    JH, D, H = get_operators_BP6(4, NY, NZ, 1.0, YN - Y0, ZN - Z0; metrics=metrics)

    print(metrics.sJ)
    # Run the simulatiom

    # Toy problem set up the final coord tranform
    dr = 2.0 / NY
    ds = 2.0 / NZ


    rc = (-1, 1, dr)
    sc = (-1, 1, ds)
    tc = (T0, TN, dt)

    run_logical(rc, sc, tc, metrics, D)
end






main()

        



