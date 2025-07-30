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
const C::Float64 = 1

struct Params{B, C, E, F}
    NR::Int64
    NS::Int64
    R_GRID::B
    S_GRID::B
    D2::C
    SAT::E
    boundary::F
end

# Exact term for Manufactured Solution - Simplified y + z + t = u(y, z, t)
function exact(t, s, r)
    # Dear god plz be merciful
    inside = C * (yf_s(r, s) + zf_r(r, s)) - t # start by removing time dependence
    return sin(inside)
end

# Exact term for Manufactured Solution Partial Derv with respect to either y or z
function exact_derv_t(t, s, r)
    inside = C * (yf_s(r, s) + zf_r(r, s)) - t
    return 0
    #return -cos(inside)
end

function exact_derv_y(t, s, r)
    inside = C * (yf_s(r, s) + zf_r(r, s)) # - t
    return C * cos(inside)
end

function exact_derv_z(t, s, r)
    inside = C * (yf_s(r, s) + zf_r(r, s)) # - t
    return 0
end

function source_term(t, s_mesh, r_mesh, res)
    ns = length(s_mesh)
    nr = length(r_mesh)
    for i in eachindex(s_mesh)
        for j in eachindex(r_mesh)
            inside = C * (yf_s(r_mesh[j], s_mesh[i]) + zf_r(r_mesh[j], s_mesh[i])) - t
            res[ns*nr + (i-1)*nr + j] += (-2*(C^2) * sin(inside))
        end
    end
    return nothing
end

function initialize(s_mesh, r_mesh, t)
    ns = length(s_mesh)
    nr = length(r_mesh)
    N = ns*nr
    res = zeros(2*N)

    for i in eachindex(s_mesh)
        for j in eachindex(r_mesh)
            res[(i-1)*nr + j] = exact(t, s_mesh[i], r_mesh[j])# U init
            res[(i-1)*nr + j + N] = exact_derv_t(t, s_mesh[i], r_mesh[j]) # v init
        end
    end
    
    return  res # ST F = U_tt - c^2Uxx
end

# Boundary Conditions:

function g_full!(R_GRID, S_GRID, NR, t, boundary)
    # Displacement U(0, Z, t)
    for i in eachindex(S_GRID)
        for j in eachindex(R_GRID)
            boundary[NR * (i - 1) + j] = exact(t, S_GRID[i], R_GRID[j])
        end
    end 
    return nothing
end

function rhs(t, x, res, ps::Params)
    # Represents RHS of ODE to solve
    # Get setup
    N = (ps.NR + 1) * (ps.NS + 1)
    
    # Step 1: Move V to U since du/dt = v
    res[1:N] = x[N+1:end]

    # Step 2: update velocity
    inter_res = zeros(N) # Tmp vector to hold boundaries for BVP
    source_term(t, ps.S_GRID, ps.R_GRID, res) # update res with source term
    SAT_Terms!(ps, x, res, t) # update res with - SAT of Boundary

    #print(res[N+1:end])
    # now solve for u everywhere
    inter_res[:] = ps.D2 \ res[N+1:end]
    res[1:N] .= inter_res
    res[1+N:end] .= inter_res
    # print("\n UGH", x)
    return nothing
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
function SAT_Terms!(ps::Params, x, res, t)
    # Dirichlet Terms
    # sat f and r 
    g_full!(ps.R_GRID, ps.S_GRID, ps.NR + 1, t, ps.boundary)
    res[((ps.NR + 1) * (ps.NS + 1))+1:end] .+= ps.SAT * ps.boundary
    return nothing
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

    # Get Coef Matrices for SAT Terms
    crr = spzeros((NR+ 1) * (NS + 1), (NR+ 1) * (NS + 1))
    diagonify(crr, metrics.crr)
    css = spzeros((NR+ 1) * (NS + 1), (NR+ 1) * (NS + 1))
    diagonify(css, metrics.css)

    @assert issparse(crr)
    @assert issparse(css)

    # Stick all the D2s together
    D2 = spzeros(NRp*NSp, NRp * NSp)
    for i in eachindex(D)
        D2 .+= D[i]
    end

    D2 .*= C^2
    D2 ./= J
    @assert issparse(D2)

    alpha_r = -13 / DR
    sat_coef1a = kron(HIs, Ir) * transpose((alpha_r .* css) + (css *  kron(D1s[1], Ir))) * Ef ./ J
    sat_coef1b = kron(HIs, Ir) * transpose((alpha_r .* css) + (css *  kron(D1s[1], Ir))) * Er ./ J

    alpha_s = -13 / DS
    sat_coef2a = kron(Is, HIr) * transpose((alpha_s .* crr) + (crr * kron(Is, D1s[2]))) * Es./ J
    sat_coef2b = kron(Is, HIr) * transpose((alpha_s .* crr) + (crr * kron(Is, D1s[2]))) * Ed ./ J

    SAT = sat_coef1a + sat_coef1b + sat_coef2a + sat_coef2b
    D2 = D2 + SAT # Account for D2 + SAT * U + boundary = 0

    @assert issparse(sat_coef1a)
    @assert issparse(sat_coef2b)

    boundary = zeros(N)
    ps = Params(NR, NS,
                R_GRID, S_GRID,  
                D2,
                SAT, boundary)

    
                 # Init for time stepping
    result = zeros(2 * (NR+ 1) * (NS + 1), NT)
    c = initialize(S_GRID, R_GRID, 0)

    g_full!(R_GRID, S_GRID, NR+1, 0, boundary)
    #print("Check C: ",  c[1:N], "\n\n")
    #print("Check B: ",  boundary .- c[1:N], "\n\n")


    @assert boundary == c[1:N]
     
    print("\nTiming for RK2:\n")
    #@time rk2_faster!(rhs, c, DT, result, T_GRID, ps)
    @time multi_step_ivbp!(rhs, c, DT, result, T_GRID, ps)
     print("\nTiming for Plotting:\n")
    @time plot_2d!(R_GRID, S_GRID, T_GRID, result, exact)

    # Only return the final time result for error
    return result[:, end]
end

function run(dy, dz, dt)
    #=
    Main is the entire simulation run, starting with physical definition then to logical

    Y and Z are physical, R and S are logical coordinates
    =#
    # Define Physical Meshes
    Y0, YN, dy = (-1, 1, dy)
    Z0, ZN, dz = (-1, 1, dz)
    T0, TN, dt = (0, 1, dt)

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

    print("\n Create Metrics: ")
    @time metrics = create_metrics_BP6(p, NZ, NY, zf_2, yf_2) # Initially do trivial one
    

    print("\n Get Ops: ")
    @time JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics, afc=false)

    print("\n")
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

    @assert issparse(Dr)
    @assert issparse(Ds)
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
converge_2D(exact; dt=1e-4, dy=0.25, dz=0.25, tc = (0, 1), yc=(-1, 1), zc = (-1, 1))

        



