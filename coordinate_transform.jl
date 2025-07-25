"""
Make Coordinate Transforms for 2D Wave
"""
# Start with a coordinate shift. Y -> R, Z -> S but squeezed by 2 and shifted down 2
function yf(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(r)
    y = zeros(i, j)

    # scaling factor
    alpha = 2
    # shift beta = 0
    beta = 0

    for row in 1:i
        for col in 1:j
            y[row, col] += alpha * (r[row, col] + beta)
        end
    end

    # Partial dervs
    yr = zeros(i, j)
    yr .+= alpha

    ys = zeros(size(s))

    return (y, yr, ys)
end

function zf(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(s)
    z = zeros(i, j)

    # scaling factor
    alpha = 4
    # shift beta = 0
    beta = 0

    for row in 1:i
        for col in 1:j
            z[row, col] += alpha * (s[row, col] + beta)
        end
    end

    # Partial dervs
    zs = zeros(i, j)
    zs .+= alpha

    zr = zeros(size(s))

    return (z, zr, zs)
end

function map_back(xf, yf, r_grid, s_grid)
    return nothing
end

##################
# Start with a coordinate shift. Y -> R, Z -> S but squeezed by 2 and shifted down 2
function zf_1(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(r)
    y = zeros(i, j)
     yr = zeros(i, j)
    

    # scaling factor
    alpha = 4.0
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

function yf_1(r,s)
    # Get size of r mesh for scaling / avg shift
    i, j = size(s)
    z = zeros(i, j)
    zs = zeros(i, j)

    # scaling factor
    alpha = 2.0
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

yf_s(r,s) = 2*s
zf_r(r,s) = 4*r