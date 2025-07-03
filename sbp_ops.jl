include("./diagonal_sbp.jl") # Updated to matching Thrase diagonal_sbp.jl

function sbp_operators(p::Int, y0::Int, yN::Int, z0::Int, zN::Int, 
    Ny::Int, Nz::Int, dy::Float64, dz::Float64)
    # Wrapper around Alex's 1D 1st Derivative operators, then makes the correct 2D for z and y in 
    # Erickson + Dunham 2014

    #Is
    Iy = sparse(Matrix{Float64}(I, Ny+1, Ny+1))
    Iz = sparse(Matrix{Float64}(I, Nz+1, Nz+1))

    y_grid = y0:dy:yN
    z_grid = z0:dz:zN

    mu = I(((Ny+1) * (Nz+1)))

    # Finalize with the E matrices
    ey0 = zeros(Ny+1)
    ey0[1] = 1

    eyn = zeros(Ny+1)
    eyn[end] = 1

    ez0 = zeros(Nz+1)
    ez0[1] = 1

    ezn = zeros(Nz+1)
    ezn[end] = 1

    Ef = kron(ey0, Iz)
    Er = kron(eyn, Iz)
    Es = kron(Iy, ez0)
    Ed = kron(Iy, ezn)

    # (D2y_test, S0_y, SN_y, HIy_test, Hy_test, r) = diagonal_sbp_D2(2, Ny; xc = (y0, yN))
    # (D2z_test, S0_z, SN_z, HIz_test, Hz_test, r) = diagonal_sbp_D2(2, Nz; xc = (z0, zN))

    # Grab Variable SBP operators from diagonal_sbp.jl
    (D2y, S0y, SNy, HIy, Hy, ry) = variable_diagonal_sbp_D2(p, Ny, 1; xc = (y0,yN))
    (D2z, S0z, SNz, HIz, Hz, rz) = variable_diagonal_sbp_D2(p, Nz, 1; xc = (z0,zN))

    # Build BS terms since Diag ^ only returns SN, S0
    BSy = SNy - S0y
    BSz = SNz - S0z

    return (D2y,  D2z, Iy, Iz, Hy, Hz, HIy, HIz, BSy, BSz, mu, Ef, Er, Es, Ed)
end


