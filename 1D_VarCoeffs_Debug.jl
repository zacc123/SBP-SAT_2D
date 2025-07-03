"""
Taking SImplified 2D Code, and adjusting for 1D on Y = 0, Z Mesh to Make sure my idea of stacking makes sense
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
    inside = zf_r(z, y)^3 + t^2
    return inside
end

# Exact term for Manufactured Solution Partial Derv with respect to either y or z
function exact_derv_y(t, y, z)
    return 3*zf_r(z, y)^2
end

function source_term(t, z_mesh)
    
    nz = length(z_mesh)
    res = zeros(nz)
  
    for z in eachindex(z_mesh)
        res[z] = 2 - (6*zf_r(z_mesh[z], 0))
    end
   
    return res
end

function initialize(z_mesh)
   
    nz = length(z_mesh)
    N = nz
    res = zeros(2*N)

    # now try u(t, y, z) = t + y
    # u(t=0) = 0 + y
    # v(t=0) = 1
   
    for j in eachindex(z_mesh)
        res[j] =  zf_r(z_mesh[j], 0)^3 # U init
        res[j + N] = 0 # v init, derv = 2t at t = 0 => 0
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
    (mu, HIr, Ir, Ef, Er,  BSr,  DR, NR,  R_GRID,
                 D2rr, D2ss, D2rs, D2sr, J, sJ, crr, css, D1s, JHi, Hr) = params
    

    # prelim setup to clean the rest up
    NRp = NR + 1
    u = x[1:NRp]
    D1r = D1s[2]
    # R Normal Terms
    alpha_r = -13 / DR

    sat_f = JHi *  transpose(( alpha_r .* crr) + (crr * D1r)) * Ef * (u[1] - g_z(R_GRID, 0,  t)[1]) 
    sat_r = JHi *  transpose(( alpha_r .* crr) + (crr * D1r)) * Er * (u[end] - g_z(R_GRID, 0,  t)[end])

    return  (sat_f .+ sat_r) 
    
end

function rhs(t, x, params)
    # Total right hand side of our ODEs
    (mu, HIr, Ir, Ef, Er,  BSr,  DR, NR,  R_GRID,
                 D2rr, D2ss, D2rs, D2sr, J, sJ, crr, css, D1s, JHi, Hr) = params
    N = (NR + 1)
    u = x[1:N]
    v = x[N+1:2*N]
    b = zeros(2 * N)
    b[1:N] = v # move u = v part
    D = D2rr
    b[N+1:2*N] = ((D ./ J) * u) .+ SAT_Terms(params, x, t) .+  source_term(t, R_GRID) # update v with sbp
    return b
end

"""
This is where everything will get a lil tricky, need to only pull out 1D operator
"""
function run_logical(rc, sc, tc, metrics, D, D1s, JHi)
    # Run the SBP SAT with Euler Post Coordinate Tranform
    # INPUTS:
        # Grid spacing info
        # rc = (r0, rN, dr)
        # sc = (s0, sN, ds)
        # tc = (t0, tN, dt)
        # metrics from the coordinate transform metrics
        # D logical 2nd order operators
        # D1s logical 1st order operators


    # Unpack and get meshes in order
    R0, RN, DR = rc
    S0, SN, DS = sc
    T0, TN, DT = tc

    R_GRID = R0:DR:RN
    S_GRID = S0:DS:SN

    NR = length(R_GRID) - 1
    NRp = NR + 1
    NS = length(S_GRID) - 1
    NSp = NS + 1
    N = NRp * NSp



    T_GRID = T0:DT:TN
    NT = length(T_GRID)

    # Get SBP Operators assuming orthonormal set up on R and S
    print("Timing for SBP OP Creation:\n")
    @time (D2s, D2r, Is, Ir, Hs, Hr, HIs, HIr, BSs, BSr, mu, Ef, Er, Es, Ed) = sbp_operators(2, S0, SN, R0, RN, NS, NR, DS, DR)
    
    # Get Coef Matrices for SAT Terms
    crr = spzeros(NRp * NSp, NRp * NSp)
    diagonify(crr, metrics.crr)
    crr = crr[1:NRp, 1:NRp]
    
    css = spzeros(NRp * NSp, NRp * NSp)
    diagonify(css, metrics.css)
    css = css[1:NRp, 1:NRp]

    # Get Jacobians and all into correct format
    J = zeros(N)
    stack!(J, metrics.J) # Stack J into same format as R, S
    J = J[1:NRp]
    # print(metrics.sJ)
    #->
    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    D2rr = D2rr[1:NRp, 1:NRp]

    (D2z, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(2, NR, 1; xc = (-2,2))
    (D2r_edge, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(2, NR, 1; xc = (-1,1))
    
    D2rr[1, :] = J .* D2r_edge[1, :]
    D2rr[end, :] = J .* D2r_edge[end, :]
    #D2rr[end, :] = J[1:3] .* ([1, -2, 1] ./ (DR^2))
    
    crri = Matrix(crr) \ I

    D2rr[end , :] =   J .* D2z[end, :] 
    D2rr[1 , :] =  J .* D2z[1, :] 

    print(D2rr[1, :])
    print("Other", J .* D2z[1, :])

    
    # Init

    result = zeros(2 * NRp, NT)
    c = initialize(R_GRID)

    er0 = zeros(NRp)
    er0[1] = 1

    ern = zeros(NRp)
    ern[end] = 1

    Ef = er0
    Er = ern

    params = (mu, HIr, Ir, Ef, Er,  BSr,  DR, NR,  R_GRID, D2rr, D2ss, D2rs, D2sr, J, metrics.sJ, crr, css, D1s, JHi, Hr) 
   
    print("\n")
    

    #@time forward_euler!(rhs, c, DT, result, T_GRID, params)
    print("\nTiming for RK2:\n")
    rk2!(rhs, c, DT, result, T_GRID, params)

    plot_1d!(R_GRID, T_GRID, result, exact)

    # Only return the final time result for error
    return result[:, end]
end

function run(dy, dz, dt)
    #=
    Main is the entire simulation run, starting with physical definition then to logical

    Y and Z are physical, R and S are logical coordinates
    =#
    # Define Physical Meshes
    Y0, YN, dy = (-1, 1, dz/2)
    Z0, ZN, dz = (-2, 2, dz)
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
    p = 2 # order of accuracy hehe

    metrics = create_metrics_BP6(p, NZ, NY, zf_1, yf_1) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics)

    JHi = JH[1:NZp, 1:NZp] ./ dz
    JHi = Matrix(JHi)\I

    # print(JHi)

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

    x = run_logical(rc, sc, tc, metrics, D, D1s, JHi)
    return x
   
end
#run(0.25, 0.25, 1e-4)
converge_1D(exact; dt=1e-4, dy=0.25, dz=0.5, tc = (0, 2), yc=(-1, 1), zc = (-1, 1))

        



