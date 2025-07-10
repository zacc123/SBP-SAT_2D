using Plots

function stack!(result_u, u_matrix)
    # Takes in the  y x z u matrix and returns a stacked vector [u1, u2, ... , un]
    y,z = size(u_matrix)
    for i in 1:y
        result_u[(i-1)*z + 1: i*z] = u_matrix[i, :]
    end
    return nothing
end

function get_const_z_vec(u::Matrix{Float64}, Ny::Int, Nz::Int, index::Int)
    # Assumes vector is stacked u = [Uy=1,z^T, Uy=2, z^T, .... , Uy=Ny+1, z^T]
    # Returns the vector len Ny for a given Z index
    result = zeros(Ny + 1)
    for i in 1:Ny+1
        result[i] = u[(i-1)*Nz + index]
    end
    return result
end

function get_r_vector!(input_uv, res_vec, NR, NS, s_index)
    # Takes in the stacked vector and returns U, V which are displacements and velocities in 2D array, Y rows, Z col
    res_vec[:] = input_uv[(s_index - 1)*(NR+1) + 1:(s_index)*(NR+1)]
    return nothing
end

function get_s_vector!(input_uv, res_vec, NR, NS, r_index)
    # Takes in the stacked vector and returns U, V which are displacements and velocities in 2D array, Y rows, Z col
    rng = r_index:NR+1:(NS + 1)*NR+r_index
    res_vec[:] = input_uv[rng]
    return nothing
end

function get_r_vectors!(input_uv, res_vec1, res_vec2, NR, NS, s_index1, s_index2)
    # Takes in the stacked vector and returns U, V which are displacements and velocities in 2D array, Y rows, Z col
    res_vec1[:] = input_uv[(s_index1 - 1)*(NR+1) + 1:(s_index1)*(NR+1)]
    res_vec2[:] = input_uv[(s_index2 - 1)*(NR+1) + 1:(s_index2)*(NR+1)]
    return nothing
end

function get_s_vectors!(input_uv, res_vec1, res_vec2, NR, NS, r_index1, r_index2)
    # Takes in the stacked vector and returns U, V which are displacements and velocities in 2D array, Y rows, Z col
    rng1 = r_index1:NR+1:(NS + 1)*NR+r_index1
    res_vec1[:] = input_uv[rng1]
    rng2 = r_index2:NR+1:(NS + 1)*NR+r_index2
    res_vec2[:] = input_uv[rng2]
    return nothing
end

function convert!(uv, y, z, u_res::Matrix{Float64}, v_res::Matrix{Float64})
    # Takes in the stacked vector and returns U, V which are displacements and velocities in 2D array, Y rows, Z cols
    NY = length(y)
    NZ = length(z)

    res_y1 = zeros(NY)
    res_y2 = zeros(NY)
    res_z1 = zeros(NZ)
    res_z2 = zeros(NZ)

    for i in eachindex(y)
        u_res[i, 1:NZ] = uv[(i - 1)*NZ + 1:(i)*NZ]
        v_res[i, 1:NZ] = uv[(i - 1)*NZ + 1 + NY*NZ:(i)*NZ + NY*NZ]

    end

    get_s_vector!(uv, res_y1, NZ-1, NY-1, 1)
    get_s_vector!(uv, res_y2, NZ-1, NY-1, NY)

    @assert res_y1 == u_res[:, 1]
    return nothing
end

function plot_with_exact(grid, num_sol, exact_func, params, labels)
    # Takes in a grid, numerical sol and plots overlayed with an exact solution
    # Pass in exact func and params for plotting there, and labels for ["Title", "X", "Y", "PNG_NAME"]
    plot(grid, num_sol, label="Numerical")
    plot!(grid, [exact_func(i, params) for i in grid], label="Exact")
    title!(labels[1])
    xlabel!(labels[2])
    ylabel!(labels[3])
    png(labels[4])
    return nothing
end

function diagonify(res_diag, matrix)
    i, j = size(matrix)
    for row in 1:i
        for col in 1:j
            indx = (i * (row - 1)) + col
            res_diag[indx, indx] = matrix[row, col]
        end
    end
    return nothing
end