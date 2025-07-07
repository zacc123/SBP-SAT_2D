#=

File meant to help with convergence testing

=#
function get_error_1d(sol, z_grid, t, exact)
    """
    Get the error of 1D grid and report the L2 Error
    Assumes Exact solution is actually 2 Vars, Y, Z but Y = 0
    """
    nz = length(z_grid)
    res = zeros(nz)

    for j in eachindex(z_grid)
        res[j] = exact(t, 0, z_grid[j])
    end
    return norm(res - sol, 2) 
end


function get_error_2d(sol, y_grid, z_grid, t, exact)
    """
    Get the error of 2D grid and report the L2 Error
    """
    ny = length(y_grid)
    nz = length(z_grid)
    res = zeros(ny * nz)

    for i in eachindex(y_grid)
        for j in eachindex(z_grid)
            res[(i-1)*nz + j] = exact(t, y_grid[i], z_grid[j])
            #print(y_grid[i], z_grid[j])
        end
    end

    return norm(res - sol, 2) 
end

function converge_2D(exact; dt=1e-4, dy=0.1, dz=0.1, tc = (0, 5), yc=(-1, -1), zc = (-1, 1))

    # Run Time and spatial convergence tests
    
    # Setup test points
    dts = [dt, dt / 2, dt / 4]
    dys = [dy, dy / 2, dy / 4, dy/8]
    dzs = [dz, dz / 2, dz / 4, dz / 8]

    # Y and Z boundaries
   

    # Start with time as sanity check
    time_errors = zeros(length(dts) + 1)
    
    # Spatial grids stay the same
    Y_GRID = yc[1]:dy:yc[2]
    Z_GRID = zc[1]:dz:zc[2]

    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    N = (NY + 1) * (NZ + 1)

    DY = dz/2
    DZ = dz
    
    # Only need to allocate once since it will keep the same dim each dt
        # Only focused on last point in time
    
    #= 
    TODO: Clean up if desired - 6/5/25
    Error is dominated by spatial for useful looks agh
    for t in eachindex(dts)
        DT = dts[t]
        print("\nRunning dt: $(DT)\n")
        res = run(DY, DZ, DT)[1:N]
        
        error = get_error_2d(res, Y_GRID, Z_GRID, tc[2], exact)
        time_errors[t] = error * sqrt(DZ*DY)
        print("\nL2 Error: $(error * sqrt(DZ*DY))\n")
        if t > 1
            print("Log2 Difference in Error: $(log2(time_errors[t-1] / time_errors[t]))")
        end
    end
    =#
    space_errors = zeros(length(dys) + 1)
    for i in eachindex(dys)
        # run with Dy = Dz = 0.125 since that works
        DY = dys[i]
        DZ = dzs[i]

        # scaling diff
        alpha_y = 1
        alpha_z = 1

        ds = DY/alpha_y
        dr = DZ/alpha_z
        
        Y_GRID = yc[1]:ds:yc[2]
        Z_GRID = zc[1]:dr:zc[2]


        NY = length(Y_GRID) - 1
        NZ = length(Z_GRID) - 1
        N = (NY + 1) * (NZ + 1)

        print("\nRunning DY: $(DY)\n")
        res = run(DY, DZ, dt)[1:N]
        error = get_error_2d(res, Y_GRID, Z_GRID, tc[2], exact)
        space_errors[i] = error * sqrt(dr*ds)
        print("\nL2 Error: $(error * sqrt(dr*ds))\n")

        if i > 1
            print("Log2 Difference in Error: $(log2(space_errors[i-1] / space_errors[i]))")
        end

    end

end

function converge_1D(exact; dt=1e-4, dy=0.1, dz=0.1, tc = (0, 5), yc=(0, 2), zc = (0, 2))
    # Setup test points
    dts = [dt, dt / 2, dt / 4]
    dzs = [dz, dz / 2, dz / 4]

    # Z boundaries
    Z_GRID = zc[1]:dz:zc[2]
    NZ = length(Z_GRID) - 1
    DZ = dz

    space_errors = zeros(length(dzs) + 1)
    for i in eachindex(dzs)
        # run with Dy = Dz = 0.125 since that works
        DZ = dzs[i]

        # scaling diff
        alpha_z = 2
        ds = DZ/alpha_z
        
        Z_GRID = zc[1]:ds:zc[2]

        NZ = length(Z_GRID) - 1
        N = NZ+1

        print("\nRunning DZ: $(DZ)\n")
        res = run(DZ, DZ, dt)[1:N]
        error = get_error_1d(res, Z_GRID, tc[2], exact)
        space_errors[i] = error * sqrt(ds)
        print("\nL2 Error: $(error * sqrt(ds))\n")

        if i > 1
            print("Log2 Difference in Error: $(log2(space_errors[i-1] / space_errors[i]))")
        end

    end

end