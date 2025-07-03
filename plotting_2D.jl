using Plots

# used to get the diff surface in the plotting
function exact_surface(sol, y_grid, z_grid, t, exact)
    ny = length(y_grid)
    nz = length(z_grid)
    res = zeros(ny * nz)
    for i in eachindex(y_grid)
        for j in eachindex(z_grid)
            res[(i-1)*nz + j] = exact(t, y_grid[i], z_grid[j])
        end
    end
    return res
end

function plot_2d!(Y_GRID, Z_GRID, T_GRID, result, exact; file_dir = "./result_imgs/")
    
    # Pass in the result and build 2 D array [y, z]
    # Plotting Grid lines in constant time
    NY = length(Y_GRID) - 1
    NZ = length(Z_GRID) - 1
    
    u = zeros(NY+1, NZ+1)
    v = zeros(NY+1, NZ+1)

    N = (NY+1) * (NZ+1)
    NT = length(T_GRID)

    convert!(result[:, 1], Y_GRID, Z_GRID, u, v)
    
    plot(Z_GRID, u[1, :], label="numerical")
    plot!(Z_GRID, [exact(T_GRID[1], Y_GRID[1], z) for z in Z_GRID],label="exact")
    png(file_dir * "Z at Y=0, T=0")

    plot(Y_GRID, u[:, 1], label="numerical")
    plot!(Y_GRID, [exact(T_GRID[1], y, Z_GRID[1]) for y in Y_GRID],label="exact")
    png("./result_imgs/Y at Z=0, T=0")

    convert!(result[:, end], Y_GRID, Z_GRID, u, v)
    plot(Z_GRID, u[1, :], label="numerical")
    plot!(Z_GRID, [exact(T_GRID[end], Y_GRID[1], z) for z in Z_GRID],label="exact")
    png("./result_imgs/Z at Y=0, T=1")

    plot(Y_GRID, u[:, 1], label="numerical")
    plot!(Y_GRID, [exact(T_GRID[end], y, Z_GRID[1]) for y in Y_GRID],label="exact")
    png("./result_imgs/Y at Z=0, T=1")

    plot(Y_GRID, u[:, end], label="numerical")
    plot!(Y_GRID, [exact(T_GRID[end], y, Z_GRID[end]) for y in Y_GRID],label="exact")
    png("./result_imgs/Y at Z=Nz, T=1")

    plot(Z_GRID, u[end, :], label="numerical")
    plot!(Z_GRID, [exact(T_GRID[end], Y_GRID[end], z) for z in Z_GRID],label="exact")
    png("./result_imgs/Z at Y=Ny, T=1")

    plot(Z_GRID, u[cld(NY+1, 2), :], label="numerical")
    plot!(Z_GRID, [exact(T_GRID[end], Y_GRID[cld(NY+1, 2)], z) for z in Z_GRID],label="exact")
    png("./result_imgs/Z at Y=0.5Ny, T=1")
    
    surface(Y_GRID, Z_GRID, result[1:N, 1])
    png("./result_imgs/Surface T 0.0")

    surface(Y_GRID, Z_GRID, result[1:N, cld(NT, 2)])
    png("./result_imgs/Surface T 0.5")

    surface(Y_GRID, Z_GRID, result[1:N, end])
    png("./result_imgs/Surface T 1.0")

    tmp = exact_surface(result[1:N, end], Y_GRID, Z_GRID, T_GRID[end], exact)

    surface(Y_GRID, Z_GRID, result[1:N, end] .- tmp)
    png("./result_imgs/DIFF Surface T 1.0")


    # Look at z = Zgrid[5], y = Ygrid[5]
    plot(T_GRID, result[4*(NZ+1) + 5, :], label="numerical")
    plot!(T_GRID, [exact(t, Y_GRID[5], Z_GRID[5]) for t in T_GRID],label="exact")
    png("./result_imgs/T at Y=0, T=0")

    return nothing
end

function plot_1d!(Z_GRID, T_GRID, result, exact; file_dir = "./result_imgs/")
    
    # Pass in the result and build 2 D array [y, z]
    # Plotting Grid lines in constant time
    NZ = length(Z_GRID) - 1
    N = (NZ+1)
    NT = length(T_GRID)
    u = result[1:N, 1]
    v = result[N+1:end, 1]

    
    plot(Z_GRID, u, label="numerical")
    plot!(Z_GRID, [exact(T_GRID[1], 0, z) for z in Z_GRID],label="exact")
    png(file_dir * "Z at Y=0, T=0")

    u = result[1:N, end]
    v = result[N+1:end, end]
    plot(Z_GRID, u, label="numerical")
    plot!(Z_GRID, [exact(T_GRID[end], 0, z) for z in Z_GRID],label="exact")
    png("./result_imgs/Z at Y=0, T=1")

    plot(T_GRID, result[1, :], label="numerical")
    plot!(T_GRID, [exact(t, 0, Z_GRID[1]) for t in T_GRID],label="exact")
    png("./result_imgs/Z=0 against T")

    return nothing
end