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

# Exact term for Manufactured Solution - Simplified y + z + t = u(y, z, t)
function exact(t, y, z)
    # Dear god plz be merciful
    inside = yf_s(z, y)^3 + zf_r(z, y)^3 + t^2
    return inside
end

# Exact term for Manufactured Solution Partial Derv with respect to either y or z
function exact_derv_y(t, y, z)
    return 3*zf_r(z, y)^2
end

function source_term(t, y_mesh, z_mesh)
    ny = length(y_mesh)
    nz = length(z_mesh)
    res = zeros(length(y_mesh) * length(z_mesh))
    for i in eachindex(y_mesh)
        for z in eachindex(z_mesh)
            res[(i-1)*nz + z] = 2 - (6*(yf_s(z_mesh[z], y_mesh[i]) + zf_r(z_mesh[z], y_mesh[i])))
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
            res[(i-1)*nz + j] = yf_s(z_mesh[j], y_mesh[i])^3 + zf_r(z_mesh[j], y_mesh[i])^3 # U init
            res[(i-1)*nz + j + N] = 0 # v init
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
    return [exact_derv_y(t, y, z) for y in y_mesh]
end

function g_z(z_mesh, y, t)
    return [exact(t, y, z) for z in z_mesh]
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
    b[N+1:2*N] = ((D ./ J) * u) .+ SAT_Terms(params, x, t) .+  source_term(t, y_mesh, z_mesh) # update v with sbp
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
    (mu, HIs, HIr, Is, Ir, Ef, Er, Es, Ed, BSs, BSr, DS, DR, NR, NS, R_GRID,
                 S_GRID, D2rr, D2ss, D2rs, D2sr, J, sJ, crr, css, D1s, JH, Hs, Hr) = params
    

    # prelim setup to clean the rest up
    NSp = NS+1
    NRp = NR + 1

    # Dirichlet Terms
    # params from E + D 2014
    beta = 1

    u_res = zeros(NSp, NRp)
    v_res = zeros(NSp, NRp)
    tmp = zeros((NSp) * (NRp) * 2)
    convert!(x, S_GRID, R_GRID, u_res, v_res)

    # fault (y=y1 case) normal is in r dir so use crr
    # R Normal Terms
    alpha_r = -13 / DR
    sat_f = (kron(HIs, Ir) * transpose((alpha_r .* css) + (css *  kron(D1s[1], Ir))) * Ef * (u_res[1, :] .- g_z(R_GRID, S_GRID[1],  t)))
    sat_r = (kron(HIs, Ir) * transpose((alpha_r .* css) + (css *  kron(D1s[1], Ir))) * Er * (u_res[end, :] .- g_z(R_GRID, S_GRID[end],  t)))

    alpha_s = -13 / DS
    sat_s = (kron(Is, HIr) * transpose((alpha_s .* crr) + (crr * kron(Is, D1s[2]))) * Es * (u_res[:, 1] .- g_y(S_GRID, R_GRID[1],  t)))
    sat_d = (kron(Is, HIr) * transpose((alpha_s .* crr) + (crr * kron(Is, D1s[2]))) * Ed * (u_res[:, end] .- g_y(S_GRID, R_GRID[end],  t)))
    
    return  (sat_f .+ sat_r .+ sat_s .+ sat_d) ./ J
    
end

function run_logical(p, rc, sc, tc, metrics, D, D1s, JH)
    # Run the SBP SAT with Euler Post Coordinate Tranform
    # INPUTS:
        # Grid spacing info
        # rc = (r0, rN, dr)
        # sc = (s0, sN, ds)
        # tc = (t0, tN, dt)
        # metrics from the coordinate transform metrics
        # D logical 2nd order operators
        # D1s logical 1st order operators

     # order for operators

    # Unpack and get meshes in order
    R0, RN, DR = rc
    S0, SN, DS = sc
    T0, TN, DT = tc

    R_GRID = R0:DR:RN
    S_GRID = S0:DS:SN

    NR = length(R_GRID) - 1
    NS = length(S_GRID) - 1
    NRp = NR + 1
    NSp = NS + 1
    N = NRp*NSp

    T_GRID = T0:DT:TN
    NT = length(T_GRID)

    # Get SBP Operators assuming orthonormal set up on R and S
    print("Timing for SBP OP Creation:\n")
    @time (D2s, D2r, Is, Ir, Hs, Hr, HIs, HIr, BSs, BSr, mu, Ef, Er, Es, Ed) = sbp_operators(p, S0, SN, R0, RN, NS, NR, DS, DR)
    
    # Get Jacobians and all into correct format
    J = zeros(N)
    stack!(J, metrics.J) # Stack J into same format as R, S
    #print(metrics.sJ)
    #->
    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    # Get Coef Matrices for SAT Terms
    crr = spzeros((NR+ 1) * (NS + 1), (NR+ 1) * (NS + 1))
    diagonify(crr, metrics.crr)
    css = spzeros((NR+ 1) * (NS + 1), (NR+ 1) * (NS + 1))
    diagonify(css, metrics.css)
    
    # Init
    result = zeros(2 * (NR+ 1) * (NS + 1), NT)
    c = initialize(S_GRID, R_GRID)

    D1r_corr =   (HIr * D1s[2]) # add in the Z_r term 
    D1s_corr =  (HIs * D1s[1])
    
    J_mat = sparse(Diagonal(J))

    edge_r = spzeros(NRp, NRp)
    edge_r[1,1] = 1
    edge_r[end, end] = 1

    derv_edge_r = spzeros(NRp, NRp)
    derv_edge_r[1,1] = -1
    derv_edge_r[end, end] = 1

    edge_s = spzeros(NSp, NSp)
    edge_s[1,1] = 1
    edge_s[end, end] = 1

    derv_edge_s = spzeros(NSp, NSp)
    derv_edge_s[1,1] = -1
    derv_edge_s[end, end] = 1
    print("\n")
    
    print("$(css[1,1])\n")
    print("$(crr[1,1])\n")
    D2rr = D2rr  +  (crr * kron(Is, derv_edge_r * D1r_corr))
    D2ss = D2ss +  (css * kron(derv_edge_s * D1s_corr , Ir))

    #D2rr = D2rr +  (crr * kron(Is, derv_edge_r * D1r_corr))
    #D2ss = D2ss  +  (css * kron(derv_edge_s * D1s_corr, Ir))

    print(Matrix(D2rr[1:NRp, 1:NRp]))
    print("\n")
    print(Matrix(D2ss[1:NRp:NRp*NSp, 1:NRp:NRp*NSp]))
    print("\n")

    # print(Matrix(D2ss))

    params = (mu, HIs, HIr, Is, Ir, Ef, Er, Es, Ed, BSs, BSr, DS, DR, NR, NS, R_GRID,
                 S_GRID, D2rr, D2ss, D2rs, D2sr, J, metrics.sJ, crr, css, D1s, JH, Hs, Hr)
   
    print("\n")
    

    #@time forward_euler!(rhs, c, DT, result, T_GRID, params)
    print("\nTiming for RK2:\n")
    rk2!(rhs, c, DT, result, T_GRID, params)

    plot_2d!(R_GRID, S_GRID, T_GRID, result, exact)

    # Only return the final time result for error
    return result[:, end]
end

function run(dy, dz, dt)
    #=
    Main is the entire simulation run, starting with physical definition then to logical

    Y and Z are physical, R and S are logical coordinates
    =#
    # Define Physical Meshes
    Y0, YN, dy = (-2, 2, dz/2)
    Z0, ZN, dz = (-4, 4, dz)
    T0, TN, dt = (0, 2, dt)

    Y_GRID = Y0:dy:YN
    Z_GRID = Z0:dz:ZN
    T_GRID = T0:dt:TN

    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    NT = length(T_GRID)

    NYp = NY + 1
    NZp = NZ + 1


    # Now make the coordinate transform
    p = 4 # order of accuracy hehe

    metrics = create_metrics_BP6(p, NZ, NY, zf_1, yf_1) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics)


    # print(JH)

    # Get 1st Derivative Operators
    Dr_tmp, _, _, _ =  diagonal_sbp_D1(p, NZ; xc = (-1, 1))
    Ds_tmp, _, _, _ =  diagonal_sbp_D1(p, NY; xc = (-1, 1))

    Ds = spzeros(NYp, NYp)
    Dr = spzeros(NZp, NZp)

    Ds[1, :] = Ds_tmp[1, :]
    Ds[end, :] = Ds_tmp[end, :]

    Dr[1, :] = Dr_tmp[1, :]
    Dr[end, :] = Dr_tmp[end, :]

    D1s = (Ds, Dr)
    
    
    # Run the simulatiom

    # Toy problem set up the final coord tranform
    ds = 2.0 / NY
    dr = 2.0 / NZ

    @assert dr == ds # leave in for now for testing


    rc = (-1, 1, dr)
    sc = (-1, 1, ds)
    tc = (T0, TN, dt)

    x = run_logical(p, rc, sc, tc, metrics, D, D1s, JH)
    return x
   
end
converge_2D(exact; dt=1e-4, dy=1.0, dz=1.0, tc = (0, 2), yc=(-1, 1), zc = (-1, 1))

        



