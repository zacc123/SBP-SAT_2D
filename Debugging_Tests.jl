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

function test_2d_operators_1D()
    """
    Test if operators match what they should from Almquist and Dunham (2021)
            under no translation in Y, and a squeeze of 2x in Z

            R is the logical dim for Z
            S is logical dim for Y 

            Matches most literature where Y/S is the stacking dim

    TEST_1: Define Z = 2r, Y = s
            - Confirm that 1D Z operator matches relation D_zz = D_rr ./ J
            - Confrim that 1D Y Operator also matches D_yy = D_ss ./ J

    TEST_2: Repeat TEST 1 but use a matrix addition + stencil to replace edges
    """

    # Setup Grids
    dy = 0.25
    dz = 0.25
    dt = 0.1

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

    NRp = NZp
    NSp = NYp
    NR = NZ
    NS = NY
    # Now make the coordinate transform
    p = 2 # order of accuracy hehe

    # Actually Build Metrics
    metrics = create_metrics_BP6(p, NZ, NY, zf_test_2x, yf_test_1x) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    # Get Operator from Existing Code
    JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics)

    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    JHi = JH[1:NZp, 1:NZp] ./ dz
    JHi = Matrix(JHi)\I

    ds = 2.0 / NY
    dr = 2.0 / NZ

    rc = (-1, 1, dr)
    sc = (-1, 1, ds)
    tc = (T0, TN, dt)

    J = zeros(NZp*NYp)
    stack!(J, metrics.J)

    # Get the phyisical operator in 1D
    rng = 1:NRp
    
    D1r, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (-1, 1)) # Get 1st Derv Operators to Get Full Comp
    D1z, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (Z0, ZN))

    (D2r_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, metrics.crr[rng]; xc = (-1, 1))
    (D2z, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, 1; xc = (Z0,ZN))

    rngs = 1 .+ NRp * (0:NS)
    D1s, HIs, Hs, sc = diagonal_sbp_D1(p, NY; xc = (-1, 1))
    D1y, HIy, Hy, sc = diagonal_sbp_D1(p, NY; xc = (Y0, YN))
    (D2s_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NY, metrics.css[rngs]; xc = (-1, 1))
    (D2y, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NY, 1; xc = (Y0,YN))
    
    # Repeat for SS but tricky to undo the K prod
    # Now check 1 D version
    D2rr_1D_Test1 = D2rr[1:NRp, 1:NRp] # this side of K prod is easy
    D2ss_1D = spzeros(NSp, NSp) # this one is tricky
    for i in 1:NSp
        D2ss_1D[i, :] = D2ss[i + (i-1)*NSp, i:NRp:NRp*NSp]
    end
    
    D1r_corr =  (J[rng] .* HIr * D1r) # add in the Z_r term 
    D1s_corr =  (J[1:NSp:NRp*NSp] .* HIs * D1s)

    D1z_corr = (HIz * D1z)
    D1y_corr = (HIy * D1y)

    # Make corrections
    D2rr_1D_Test1[1, :] = D1r_corr[1, :]
    D2rr_1D_Test1[end, :] = -1 .* D1r_corr[end, :]
    
    D2z[1, :] = D1z_corr[1, :]
    D2z[end, :] = -1 .* D1z_corr[end, :]

    D2y[1, :] = D1y_corr[1, :]
    D2y[end, :] = -1 .* D1y_corr[end, :]
    
    D2ss_1D_Test1 = deepcopy(D2ss_1D)
    D2ss_1D_Test1[1, :] = D1s_corr[1, :]
    D2ss_1D_Test1[end, :] = -1 .* D1s_corr[end, :]
    
    # TEST_1:
    # Now check if relation is true
    @assert D2z == D2rr_1D_Test1 ./ J[1:NZp]

    #print(Matrix(D2y))
    #print("\n\n")
    #print(Matrix(D2ss_1D_Test1 ./ J[1:NSp:NRp*NSp]))
    @assert D2y == D2ss_1D_Test1 ./ J[1:NSp:NRp*NSp]

    # TEST_2 See if we can do the same as 1 but with some matrix addition instead of 
        # indexing
    Ir = sparse(Diagonal(ones(NRp)))
    Is = sparse(Diagonal(ones(NSp)))

    edge_r = spzeros(NRp, NRp)
    edge_r[1,1] = 1
    edge_r[end, end] = 1

    derv_edge_r = spzeros(NRp, NRp)
    derv_edge_r[1,1] = 1
    derv_edge_r[end, end] = -1

    edge_s = spzeros(NSp, NSp)
    edge_s[1,1] = 1
    edge_s[end, end] = 1

    derv_edge_s = spzeros(NSp, NSp)
    derv_edge_s[1,1] = 1
    derv_edge_s[end, end] = -1

    D2rr_1D_Test2 = D2rr[1:NRp, 1:NRp]
    D2rr_1D_Test2 = D2rr_1D_Test2 - (edge_r * D2rr_1D_Test2) + (derv_edge_r * D1r_corr)

    D2ss_1D_Test2 = deepcopy(D2ss_1D)
    D2ss_1D_Test2 = D2ss_1D_Test2 - (edge_s * D2ss_1D_Test2) + (derv_edge_s * D1s_corr)
    # print(Array(D1rr_1D_Test2))
    # print(Array(D2rr_1D))

    # TEST_2
    @assert D2rr_1D_Test2 == D2rr_1D_Test1
    @assert D2ss_1D_Test2 == D2ss_1D_Test1

    return D2z, D2rr
end

function test_2d_operators_2D()
    """
    Test if operators match what they should from Almquist and Dunham (2021)
            under no translation in Y, and a squeeze of 2x in Z

            R is the logical dim for Z
            S is logical dim for Y 

            Matches most literature where Y/S is the stacking dim

            We'll Now do stuff to copy 2D operators

    TEST_1: Define Z = r, Y = s
            - Check that basic operators match the Erickson + Dunham relations
                kron(Is, D2r) = Drr, 
                kron(D2s, Ir) = Dss, 

    TEST_2: Define Z = 2r, Y = s
            - Confirm that 1D Z operator matches relation D_ = D_rr ./ J
            - Confrim that 1D Y Operator also matches D_yy = D_ss ./ J

    """

    # Setup Grids
    dy = 0.25
    dz = 0.25
    dt = 0.1

    Y0, YN, dy = (-1, 1, dz)
    Z0, ZN, dz = (-1, 1, dz)
    T0, TN, dt = (0, 2, dt)

    Y_GRID = Y0:dy:YN
    Z_GRID = Z0:dz:ZN
    T_GRID = T0:dt:TN

    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    NT = length(T_GRID)

    NYp = NY + 1
    NZp = NZ + 1

    NRp = NZp
    NSp = NYp
    NR = NZ
    NS = NY
    # Now make the coordinate transform
    p = 2 # order of accuracy hehe

    # Actually Build Metrics
    metrics = create_metrics_BP6(p, NZ, NY, zf_test_1x, yf_test_1x) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    # Get Operator from Existing Code
    JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics)

    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    ds = 2.0 / NY
    dr = 2.0 / NZ

    rc = (-1, 1, dr)
    sc = (-1, 1, ds)
    tc = (T0, TN, dt)

    J = zeros(NZp*NYp)
    stack!(J, metrics.J)

    # Get the phyisical operator in 1D

    rng = 1:NRp

    D1r, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (-1, 1)) # Get 1st Derv Operators to Get Full Comp
    D1z, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (Z0, ZN))

    (D2r_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, metrics.crr[rng]; xc = (-1, 1))
    (D2z, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, 1; xc = (Z0,ZN))

    rngs = 1 .+ NRp * (0:NS)

    D1s, HIs, Hs, sc = diagonal_sbp_D1(p, NY; xc = (-1, 1))
    D1y, HIy, Hy, sc = diagonal_sbp_D1(p, NY; xc = (Y0, YN))
    (D2s_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NY, metrics.css[rngs]; xc = (-1, 1))
    (D2y, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NY, 1; xc = (Y0,YN))
    
    # Repeat for SS but tricky to undo the K prod
    # Now check 1 D version
    D2rr_1D_Test1 = D2rr[1:NRp, 1:NRp] # this side of K prod is easy
    D2ss_1D = spzeros(NSp, NSp) # this one is tricky
    for i in 1:NSp
        D2ss_1D[i, :] = D2ss[i + (i-1)*NSp, i:NRp:NRp*NSp]
    end

    D1r_corr =  (J[rng] .* HIr * D1r) # add in the Z_r term 
    D1s_corr =  (J[1:NSp:NRp*NSp] .* HIs * D1s)

    D1z_corr = (HIz * D1z)
    D1y_corr = (HIy * D1y)

    # Make corrections
    D2rr_1D_Test1[1, :] = D1r_corr[1, :]
    D2rr_1D_Test1[end, :] = -1 .* D1r_corr[end, :]
    
    D2z[1, :] = D1z_corr[1, :]
    D2z[end, :] = -1 .* D1z_corr[end, :]

    D2y[1, :] = D1y_corr[1, :]
    D2y[end, :] = -1 .* D1y_corr[end, :]
    
    D2ss_1D_Test1 = deepcopy(D2ss_1D)
    D2ss_1D_Test1[1, :] = D1s_corr[1, :]
    D2ss_1D_Test1[end, :] = -1 .* D1s_corr[end, :]
    
    # TEST_0:
    # Now check if relation is true
    @assert D2z == D2rr_1D_Test1 ./ J[1:NZp]
    @assert D2y == D2ss_1D_Test1 ./ J[1:NSp:NRp*NSp]

    # TEST_2 See if we can do the same as 1 but with some matrix addition instead of 
        # indexing
    Ir = sparse(Diagonal(ones(NRp)))
    Is = sparse(Diagonal(ones(NSp)))

    edge_r = spzeros(NRp, NRp)
    edge_r[1,1] = 1
    edge_r[end, end] = 1

    derv_edge_r = spzeros(NRp, NRp)
    derv_edge_r[1,1] = 1
    derv_edge_r[end, end] = -1

    edge_s = spzeros(NSp, NSp)
    edge_s[1,1] = 1
    edge_s[end, end] = 1

    derv_edge_s = spzeros(NSp, NSp)
    derv_edge_s[1,1] = 1
    derv_edge_s[end, end] = -1

    D2rr_1D_Test2 = D2rr[1:NRp, 1:NRp]
    D2rr_1D_Test2 = D2rr_1D_Test2 - (edge_r * D2rr_1D_Test2) + (derv_edge_r * D1r_corr)

    D2ss_1D_Test2 = deepcopy(D2ss_1D)
    D2ss_1D_Test2 = D2ss_1D_Test2 - (edge_s * D2ss_1D_Test2) + (derv_edge_s * D1s_corr)
    # print(Array(D1rr_1D_Test2))
    # print(Array(D2rr_1D))

    # TEST_0
    @assert D2rr_1D_Test2 == D2rr_1D_Test1
    @assert D2ss_1D_Test2 == D2ss_1D_Test1

    # TEST_1 
    # No longer works because we're adjusting the Boundary terms on the physical as well
    @assert kron(Is, D2rr_1D_Test1) == kron(Is, D2z) ./ J
    @assert kron(D2ss_1D_Test1, Ir) == kron(D2y, Ir) ./ J
    # @assert kron(D2ss_1D_Test2, Ir) == kron(D2y, Ir) ./ J

    # TEST_2 Add coordinate transform

    # Setup Grids
    dy = 0.125
    dz = 0.25
    dt = 0.1

    Y0, YN, dy = (-1, 1, dy)
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

    NRp = NZp
    NSp = NYp
    NR = NZ
    NS = NY
    # Now make the coordinate transform
    p = 2 # order of accuracy hehe

    # Actually Build Metrics
    metrics = create_metrics_BP6(p, NZ, NY, zf_test_2x, yf_test_1x) # Initially do trivial one
    #metrics = create_metrics_BP6(4, NY, NZ) # Initially do trivial one

    # Get Operator from Existing Code
    JH, D, H = get_operators_BP6(p, NZ, NY, 1.0, ZN - Z0, YN - Y0; metrics=metrics)

    D2rr, D2ss, D2rs, D2sr = D # Unpack D2 ops plus cross terms : , /

    ds = 2.0 / NY
    dr = 2.0 / NZ

    rc = (-1, 1, dr)
    sc = (-1, 1, ds)
    tc = (T0, TN, dt)

    J = zeros(NZp*NYp)
    stack!(J, metrics.J)

    D1r, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (-1, 1)) # Get 1st Derv Operators to Get Full Comp
    D1z, HIr, Hr, rc = diagonal_sbp_D1(p, NZ; xc = (Z0, ZN))

    # Get the phyisical operator in 1D
    rng = 1:NRp
    (D2r_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, metrics.crr[rng]; xc = (-1, 1))
    (D2z, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, 1; xc = (Z0,ZN))

    D1s, HIs, Hs, sc = diagonal_sbp_D1(p, NY; xc = (-1, 1))
    D1y, HIy, Hy, sc = diagonal_sbp_D1(p, NY; xc = (Y0, YN))
    rngs = 1 .+ NRp * (0:NS)
    (D2s_1d, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, metrics.css[rngs]; xc = (-1, 1))
    (D2y, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, NZ, 1; xc = (Y0,YN))
    
    # Repeat for SS but tricky to undo the K prod
    # Now check 1 D version
    D2rr_1D_Test1 = D2rr[1:NRp, 1:NRp] # this side of K prod is easy
    D2ss_1D = spzeros(NSp, NSp) # this one is tricky
    for i in 1:NSp
        D2ss_1D[i, :] = D2ss[i + (i-1)*NSp, i:NRp:NRp*NSp]
    end
    
    D1r_corr =  (J[rng] .* HIr * D1r) # add in the Z_r term 
    D1s_corr =  (J[1:NSp:NRp*NSp] .* HIs * D1s)

    D1z_corr = (HIz * D1z)
    D1y_corr = (HIy * D1y)

    D2z[1, :] = D1z_corr[1, :]
    D2z[end, :] = -1 .* D1z_corr[end, :]

    D2y[1, :] = D1y_corr[1, :]
    D2y[end, :] = -1 .* D1y_corr[end, :]

    # TEST_2 See if we can do the same as 1 but with some matrix addition instead of 
        # indexing
    Ir = sparse(Diagonal(ones(NRp)))
    Is = sparse(Diagonal(ones(NSp)))

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

    D2rr_1D_Test2 = D2rr[1:NRp, 1:NRp]
    D2rr_1D_Test2 = D2rr_1D_Test2 - (edge_r * D2rr_1D_Test2) + (derv_edge_r * D2r_1d)

    D2ss_1D_Test2 = deepcopy(D2ss_1D)
    D2ss_1D_Test2 = D2ss_1D_Test2 - (edge_s * D2ss_1D_Test2) + (derv_edge_s * D2s_1d)
    # print(Array(D1rr_1D_Test2))
    # print(Array(D2rr_1D))

    # TEST_1
    # Remove the actual edges and replace
    D2rr = D2rr - (kron(Is, edge_r) * D2rr) +  kron(Is, derv_edge_r * D1r_corr)
    D2ss = D2ss - (kron(edge_s, Ir) * D2ss) +  kron(derv_edge_s * D1s_corr, Ir)

    #print(Array(D2rr))
    #print("\n")
    #print(Array(kron(Is, D2z) .* J))
    @assert D2rr == kron(Is, D2z) .* J
    @assert D2ss == kron(D2y, Ir) .* J
    
    # @assert kron(D2ss_1D_Test2, Ir) == kron(D2y, Ir) ./ J
    return D2z, D2rr
end

#### Extras to run in tests
function zf_test_2x(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(r)
    y = zeros(i, j)
     yr = zeros(i, j)
    

    # scaling factor
    alpha = 2.0
    # shift beta = 0
    beta = 0.0
    charlie = 2

    for row in 1:i
        for col in 1:j
            y[row, col] +=  alpha * r[row, col]
            yr[row, col]+= alpha
        end
    end

    # Partial dervs
    ys = zeros(size(s))

    return (y, yr, ys)
end

function zf_test_1x(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(r)
    y = zeros(i, j)
     yr = zeros(i, j)
    

    # scaling factor
    alpha = 1.0
    # shift beta = 0
    beta = 0.0
    charlie = 2

    for row in 1:i
        for col in 1:j
            y[row, col] +=  alpha * r[row, col]
            yr[row, col]+= alpha
        end
    end

    # Partial dervs
    ys = zeros(size(s))

    return (y, yr, ys)
end

function yf_test_1x(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(s)
    z = zeros(i, j)
    zs = zeros(i, j)

    # scaling factor
    alpha = 1.0
    # shift beta = 0
    beta = 0.0
    charlie = 2

    for row in 1:i
        for col in 1:j
            z[row, col] += alpha * s[row, col]
            zs[row, col] += alpha
        end
    end

    # Partial dervs
    
    

    zr = zeros(size(s))

    return (z, zr, zs)
end

yf_s(r,s) = 1*s
zf_r(r,s) = 2*r


function main()
    test_2d_operators_1D()
    test_2d_operators_2D()
end

main()