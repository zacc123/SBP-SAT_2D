# Define a forward euler function
function forward_euler!(f, c, h, result, time_mesh, params)
    # Forward Euler function to find an approximate solution for an ODE
    # INPUTS
        # f is a function fo t and y, f(t, y), for the ODE, y is a vector
        # c is an initial condition vector
            # -> We'll assume its of dim 2(N+1) since we have to keep track of U_1 and U_2
        # h is a step size in time
        # result is our resulting matrix of dim (N+1xdim(y))
        # time_mesh is a 1D array for our time scale with N+1 nodes
        # params is anything needed in f
    # Output: a resultant vector with the solution

    # Start with initial condition
    result[:, 1] = c
    for n = 2:length(time_mesh)
        result[:, n] = result[:, n-1] + h*f(time_mesh[n-1], result[:, n-1], params)
    end
    
    return nothing
end

function rk2!(f, c, h, result, time_mesh, params)
    # Explicit Trapezoid function to find an approximate solution for an ODE
    # INPUTS
        # f is a function fo t and y, f(t, y), for the ODE, y is a vector
        # c is an initial condition vector
            # -> We'll assume its of dim 2(N+1) since we have to keep track of U_1 and U_2
        # h is a step size in time
        # result is our resulting matrix of dim (N+1xdim(y))
        # time_mesh is a 1D array for our time scale with N+1 nodes
        # params is anything needed in f
    # Output: a resultant vector with the solution

    # Start with initial condition
    result[:, 1] = c
    for n = 2:length(time_mesh)
        result[:, n] = result[:, n-1] + h*f(time_mesh[n-1] + 0.5*h, result[:, n-1] + 0.5*h*f(time_mesh[n-1], result[:, n-1], params), params)
    end
    
    return nothing
end
#
function rk2_faster!(f!, c, h, result, time_mesh, params)
    # Explicit Trapezoid function to find an approximate solution for an ODE
    # INPUTS
        # f is a function fo t and y, f(t, y), for the ODE, y is a vector
        # c is an initial condition vector
            # -> We'll assume its of dim 2(N+1) since we have to keep track of U_1 and U_2
        # h is a step size in time
        # result is our resulting matrix of dim (N+1xdim(y))
        # time_mesh is a 1D array for our time scale with N+1 nodes
        # params is anything needed in f
    # Output: a resultant vector with the solution

    # Start with initial condition
    result[:, 1] = c

    # Goal is to not have to allocate anymore memory 
    tmp_vec1 = zeros(length(c))
    tmp_vec2 = zeros(length(c))

    for n = 2:length(time_mesh)
        # Do the 2 Step RK

        # Step1
        f!(time_mesh[n-1], result[:, n-1], tmp_vec1, params)
        tmp_vec1 = ((0.5*h) .* tmp_vec1)  .+ result[:, n-1]

        # Step2
        f!(time_mesh[n-1] + 0.5*h, tmp_vec1, tmp_vec2, params)
        result[:, n] = result[:, n-1] .+ (h .* tmp_vec2)
    end
    
    return nothing
end
