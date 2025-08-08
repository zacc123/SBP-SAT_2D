using LinearAlgebra
using SparseArrays
using NaNMath

function sbp_operators_2D(p::Int, y0::Int, yN::Int, z0::Int, zN::Int, 
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

    Ef = kron(ey0*ey0', Iz)
    Er = kron(eyn*eyn', Iz)
    Es = kron(Iy, ez0*ez0')
    Ed = kron(Iy, ezn*ezn')

    @assert issparse(Ef)
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



function sbp_operators_1D(SBPp::Int, x0::Float64, xN::Float64, Nx::Int, dx::Float64, μ::Float64, metrics)
    # Wrapper around Alex's 1D 1st Derivative operators, then makes the correct 2D for z and y in 
    # Erickson + Dunham 2014

    AFC = true # Use Adapted Fully Compatible Operators from A&D

    # Set up Initial Bits
    Ix = sparse(Matrix{Float64}(I, Nx+1, Nx+1))

    x_grid = x0:dx:xN
    Lx = xN - x0

    # Finalize with the E matrices
    ex0 = spzeros(Nx+1)
    ex0[1] = 1

    exn = spzeros(Nx+1)
    exn[end] = 1

    Es = ex0*ex0'
    Ed = exn*exn'

    @assert issparse(Es)
 
    # Grab Variable SBP operators from diagonal_sbp.jl
    print("\n Get D2 Ops Almquist: ")
    @time JH, D, H = get_operators_BP6(SBPp, Nx, Nx, μ, Lx, Lx; metrics=metrics, afc=AFC)
    # Get 1st Derv Operators for SAT Terms and Traction
    
    print("\n")
    D1x, HIx, Hx, _ =  diagonal_sbp_D1(SBPp, Nx; xc = (-1, 1))

    D1x_SAT = spzeros(Nx+1, Nx+1) # Set normal for traction
    D1x_SAT[1, :] = D1x[1, :]
    D1x_SAT[end, :] = D1x[end, :]

    T = [μ .* D1x] # Traction term

    # Now get ready for SAT

    # Following shennanigans to account for get_ops in 2D, but I want this in 1D
    # Grab the jacobian for only first R Row of S, making it a 1D problem
    J_tmp = zeros((Nx+1)*(Nx+1))
    J = zeros(Nx+1)
    stack!(J_tmp, metrics.J) # Stack J into same format as R, S
    J .= J_tmp[1:Nx+1]

    # Get Coef Matrices for SAT Terms
    crr_tmp = spzeros((Nx+1) * (Nx+1),(Nx+1) * (Nx+1))
    diagonify(crr_tmp, metrics.crr)
    crr = crr_tmp[1:Nx+1, 1:Nx+1]
    @assert issparse(crr)

    D2 = D[1][1:Nx+1, 1:Nx+1]
    D2 ./= J

    @assert issparse(D2)

    # Build Coefficients ahead of time since these will be multiplied by U - Boundary
    alpha_r = -13 / dx
    sat_coef1a = HIx * transpose((alpha_r .* crr) + (crr *  D1x_SAT)) * Es ./ J
    sat_coef1b = HIx * transpose((alpha_r .* crr) + (crr *  D1x_SAT)) * Ed ./ J

    #print(Array(sat_coef1a))
    # Account for QD case where 0 = D2u + SATu + f
    SAT = sat_coef1a + sat_coef1b
    D2 += SAT

    #if AFC == false || SBPp > 2
        # Issue with LU for SBPp == 2 AFC because first and last rows are all 0s
         D2 = lu(D2)
    #end

    B = [sat_coef1a, sat_coef1b]
    e = (Es, Ed)

    return (D2 = D2,
            B = B,
            T = T,
            e = e)
end
# Changed to only include 2 Faces
function bdry_vec_strip!(g, B, x, slip_data, remote_data, Lx)

    g[:] .= 0

    # fault (Dirichlet):
    vf = slip_data
    g[:] += B[1] * vf


    # FACE 2 (Dirichlet):
    vf = remote_data
    g[:] += B[2] * vf

    return nothing
  
end

function computetraction_stripped(T, u, e)
    e1 = e[1]
    return (e1' * T[1]* u)[1] 
end

function rateandstate(V, psi, σn, ϕ, η, a, V0)
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)
    dfdV  = a .* (1 ./ sqrt.(1 .+ (V .* Y).^2)) .* Y
  
    g    = σn .* f    .+ η .* V .- ϕ
    dgdV = σn .* dfdV .+ η
    #print(g, dgdV)
    (g[1], dgdV[1])
end

function newtbndv(func, xL, xR, x; ftol = 1e-6, maxiter = 500, minchange=0,
                    atolx = 1e-4, rtolx = 1e-4)
    (fL, _) = func(xL)
    (fR, _) = func(xR)
    if fL .* fR > 0
        #print(fL, "\n",fR, "\n")
      return (typeof(x)(NaN), typeof(x)(NaN), -maxiter)
    end
  
    (f, df) = func(x)
    dxlr = xR - xL
  
    for iter = 1:maxiter
      dx = -f / df
      x  = x + dx
  
      if x < xL || x > xR || abs(dx) / dxlr < minchange
        x = (xR + xL) / 2
        dx = (xR - xL) / 2
      end
  
      (f, df) = func(x)
  
      if f * fL > 0
        (fL, xL) = (f, x)
      else
        (fR, xR) = (f, x)
      end
      dxlr = xR - xL
  
      if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x))
        return (x, f, iter)
      end
    end
    return (x, f, -maxiter)
  end


  ### UTILS
using Plots
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using DifferentialEquations
using Interpolations


function interp1(xpt, ypt, x)

  knots = (xpt,) 
  itp = interpolate(knots, ypt, Gridded(Linear()))
  #itp[x]  # endpoints of x must be between xpt[1] and xpt[end]
end
      
function create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, t, RSVinit, δ, τz0, θ)


  path_to_slip = pth * "slip.dat"
  # slip.dat is a file that stores time, max(V) and slip at all the stations:
  open(path_to_slip, "w") do io
    write(io,"0.0 0.0 ")
      for i in 1:length(flt_loc)
        write(io,"$(flt_loc[i]) ")
      end
        write(io,"\n")
    end
  
  #write out initial data into devol.txt:
  vv = Array{Float64}(undef, 1, 2+length(flt_loc))
    vv[1] = t
    vv[2] = NaNMath.log10(RSVinit)
    vv[3:end] = δ[flt_loc_indices]
    open(path_to_slip, "a") do io
        writedlm(io, vv)
    end

  # write out initial data into station files:

  # fltst_dpXXX.txt is a file that stores time and time-series of slip, log10(slip_rate), 
  # shear_stress and log10(state) at depth of z = XXX km, where XXX is each of the fault station depths.
  # First we write out initial data into each fltst_dpXXX.txt:

  for n = 1:length(station_strings)
    XXX = pth * "fltst_strk"*station_strings[n]*".txt"
    ww = Array{Float64}(undef, 1, 5)
    ww[1] = t
    ww[2] = δ[station_indices[n]]
    ww[3] = NaNMath.log10(RSVinit)
    ww[4] = τz0
    ww[5] = NaNMath.log10(θ[station_indices[n]])  # 
    open(XXX, "w") do io
      write(io, "# problem=SEAS Benchmark BP1-QD\n")  # 
      write(io, "# code=Thrase\n")
      write(io, "# modeler=B. A. Erickson\n")
      write(io, "# date=2023/01/09\n")
      write(io, "# element size=xx m\n")
      write(io, "# location=on fault, z = "*string(parse(Int64, station_strings[n])/10)*" km\n")
      write(io, "# Lz = 40 km\n")
      write(io, "t slip slip_rate shear_stress state\n")

      writedlm(io, ww)
    end
  end

end

function write_to_file(pth, ψδ, t, i, zf, flt_loc, flt_loc_indices, station_strings, station_indices, p, base_name="", tdump=100)
  
  path_to_slip = pth * "slip.dat"
  Vmax = 0.0

  if isdefined(i,:fsallast) 
    δNp = 1
    Nz = 0
    dψV = i.fsallast
    dψ = @view dψV[1:δNp]
    V = @view dψV[δNp .+ (1:Nz+1)]
    Vmax = maximum(abs.(extrema(V)))
    δ = @view ψδ[δNp .+ (1:Nz+1)]
    ψ = @view ψδ[1:δNp]
    τf = p.τf
  
 
    θ = (p.RSDc * exp.((ψ .- p.RSf0) ./ p.RSb)) / p.RSV0  # Invert ψ for θ.
  
    if mod(ctr[], p.save_stride_fields) == 0 || t == (p.sim_years ./ 31556926)
      vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = NaNMath.log10(Vmax)
      vv[3:end] = δ[flt_loc_indices]
      open(path_to_slip, "a") do io
        writedlm(io, vv)
      end

      for i = 1:length(station_indices)
        ww = Array{Float64}(undef, 1, 5)
        ww[1] = t
        ww[2] = δ[station_indices[i]]
        ww[3] = NaNMath.log10(V[station_indices[i]])
        ww[4] = τf[station_indices[i]]
        ww[5] = NaNMath.log10(θ[station_indices[i]])

        XXX = pth * "fltst_strk"*station_strings[i]*".txt"
        open(XXX, "a") do io
            writedlm(io, ww)
        end
      end
      
    end
  
    global ctr[] += 1
  end

  Vmax
end




# find_ind() differentiates b/t phases by defining
# interseismic when max slip rate < 10^-3 m/s
# mv is maximum slip rate (log10 m/s) 
function find_ind(mv)
  ind = [1]
  int = 1
  cos = 0
  for i = 2:length(mv)
    if mv[i] > -3 && int == 1 && cos == 0
      append!(ind, i);
      int = 0;
      cos = 1;
    end
  
    if mv[i] < -3 && int == 0 && cos == 1
      append!(ind, i-1)
      int = 1
      cos = 0
    end
  end


  ind = append!(ind, length(mv));  #tack on for plotting any part of an incomplete coseismic/interseismic phase
  
  return ind
end

# plot_slip will plot slip contours from devol.txt - every 5 years in blue during interseismic, 
# every 1 second in red during coseismic
function plot_slip(filename)

  grid = readdlm(filename, Float64)
  sz = size(grid)
  print("SIZE: $(sz)")
  flt_loc = grid[1,3:end]
  T = grid[2:sz[1],1]
  maxV = grid[2:end, 2]
  slip = grid[2:sz[1], 3:3]
  N = size(slip)[2]


  ind = find_ind(maxV);        #finds indices for inter/co-seismic phases
  interval = [5*31556926 1]   #plot every 5 years and every 1 second
  
  ct = 0   #this counts the number of events


  #Assumes an initial interseismic period
  #This for-loop only plots completed phases
  print("IND: $(ind)")
  for i = 1:2:length(ind)-2
    
    T1 = T[ind[i]]:interval[1]:T[ind[i+1]];
    print("\nT1: $(T1)\n")
    W1 = interp1(T,slip[:,1],T1)';
    
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    if i == 1
      plot(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    else
      plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    end

   
    T1 = T[ind[i+1]]:interval[2]:T[ind[i+2]];


    W1 = interp1(T,slip[:,1],T1)';
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    plot!(W1, -flt_loc, linecolor = :red, legend = false) #interseismic phase

    ct = ct+1;
  end

  
  # plot remainder of an incomplete interseismic period:
  i = length(ind)-1;
  T1 = T[ind[i]]:interval[1]:T[ind[i+1]];
  W1 = interp1(T,slip[:,1],T1)';
      print("\nT1 2: $(T1)\n")
      print("\nW1 2: $(W1)\n")
      print("\ni: $(i)\n")
      print("\nflt loc 2: $(flt_loc)\n")
      nodes = length(W1)
      for j = 2:N 
        w1 = interp1(T,slip[:,j],T1)';
        W1 = [W1; w1]
      end
      print("\nW1 3: $(W1)\n")
      if i == 1
        #plot(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
        plot(1:nodes, W1, linecolor = :blue, legend = false) #interseismic phase
      else
        plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
      end

      xlabel!("Cumulative Slip (m)")
      ylabel!("Depth (km)")
end


function plot_global(field, filename)

  @show filename
  grid = readdlm(filename)#, Float64)  # some elements cannot be parsed as numbers, 
                                       # a heterogeneous array of numbers and strings is returned.
  sz = size(grid)
  
  # indexing `grid` starting at row 9 to skip the header info
  T = grid[11:end, 1]  # Get time.
  T = T ./ 31556926 # convert to years.
 @show field
  if field == "maxV"
    y = grid[11:sz[1],2]
    plot(T, y)
  elseif field == "moment_rate"
    y = grid[11:sz[1],3]
    plot(T, y)
  else
    print("field not recognized")
  end
  gui()
    #return nothing
end



# plot_fault_time_series will plot field "field" from "filename".
# "field" has to be one of "slip", "V", "shear_stress", "state"
function plot_fault_time_series(field, filename)

  @show filename
  grid = readdlm(filename)#, Float64)  # some elements cannot be parsed as numbers, 
                                       # a heterogeneous array of numbers and strings is returned.
  sz = size(grid)
  
  # indexing `grid` starting at row 9 to skip the header info
  T = grid[9:end, 1]  # Get time.
  T = T ./ 31556926 # convert to years.
 @show field
  if field == "slip"
    y = grid[9:sz[1],2]
    plot(T, y, label="slip")
    ylabel!("slip [m]")
  elseif field == "slip_rate"
    y = grid[9:sz[1],3]
    plot(T, y, label="slip rate")
    ylabel!("slip rate [m/s]")
  elseif field == "shear_stress"
    y = grid[9:sz[1],4]
    plot(T, y)
    ylabel!("shear stress [MPa]")
  elseif field == "state"
    y = grid[9:sz[1],7]
    plot(T, y)
    ylabel!("state")
  else
    print("field not recognized")
  end
  xlabel!("time [yr]")
  gui()
    #return nothing
end

# Function for reading in numerical parameters 
function read_params(f_name)
  f = open(f_name, "r")
  tmp_params = []
  while ! eof(f)
      s = readline(f)
      if s[1] != '#'
          push!(tmp_params, split(s, '=')[2])
          flush(stdout)
      end
  end
  close(f)

  #(pth, stride_space, stride_time, xc, zc, Nx, Nz,
  #sim_years, Vp, ρ, cs, σn, RSamin, RSamax, RSb, RSDc,
  #RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])


    params = Vector{Any}(undef, 23)
    params[1] = strip(tmp_params[1])
    
    params[2] = parse(Int64, tmp_params[2])
    params[3] = parse(Int64, tmp_params[3])
    params[4] = (parse(Float64, tmp_params[4]), parse(Float64, tmp_params[5]))
    params[5] = (parse(Float64, tmp_params[6]), parse(Float64, tmp_params[7]))
    params[6] = parse(Int64, tmp_params[8])
    params[7] = parse(Int64, tmp_params[9])
    for i = 10:length(tmp_params)-1
      params[i-2] = parse(Float64, tmp_params[i])
    end

    params[23] = parse(Int64, tmp_params[25])

  return params
end


function read_params_BP6(f_name)
  f = open(f_name, "r")
  tmp_params = []
  while ! eof(f)
      s = readline(f)
      if s[1] != '#'
          push!(tmp_params, split(s, '=')[2])
          flush(stdout)
      end
  end
  close(f)

 
  #(pth, stride_space, stride_time, xc, zc, Hx, Hz, Nr, Ns, dx, dz, el_r, el_s,
    #sim_years, Vp, ρ, cs, σn_0, RSa, RSb, RSD_RS,
    #RSf0, RSV0, RSVinit, RSLf, lz, μshear, τ_init, η, q_0, 
    #t_off, α, β, φ, k, η_visc, state_law, SBPp) = read_params(localARGS[1])

    params = Vector{Any}(undef, 35)
    params[1] = strip(tmp_params[1])
    
    params[2] = parse(Int64, tmp_params[2])
    params[3] = parse(Int64, tmp_params[3])
    params[4] = (parse(Float64, tmp_params[4]), parse(Float64, tmp_params[5]))
    params[5] = (parse(Float64, tmp_params[6]), parse(Float64, tmp_params[7]))
    params[6] = parse(Int64, tmp_params[8])
    params[7] = parse(Int64, tmp_params[9])
    params[8] = parse(Int64, tmp_params[10])
    params[9] = parse(Int64, tmp_params[11])
    for i = 12:length(tmp_params)-2
      params[i-2] = parse(Float64, tmp_params[i])
    end
    params[34] = tmp_params[36]

    params[35] = parse(Int64, tmp_params[37])

  return params
end

# animate_slip will plot slip profiles against depth for every time step computed:
function animate_slip(S, δNp, zf, stride_time)

  m = length(zf)
  no_time_steps = size(S.t)
  slip_final = S.u[end][end]

  for i = 1:stride_time:no_time_steps[1]

    slip_t = S.u[i][δNp+1:end] # slip at time t
    #pyplot()
    display(plot(slip_t, -zf, xtickfont=font(18),
    ytickfont=font(18),
    guidefont=font(18),
    legendfont=font(18), ylabel = "Depth (km)", xlabel = "Slip (m)", xlims = (0, slip_final)))
    sleep(0.1)
  end

  #nothing
end







function write_to_file_BP6(pth, ψδ, t, i, zf,flt_loc, flt_loc_indices, stations, station_indices, p, μshear, dz, base_name="", tdump=100)
  
  path_to_global = pth * "global.dat"
  path_to_slip = pth * "slip.dat"

  Vmax = 0.0

  if isdefined(i,:fsallast) 
    δNp = p.δNp
    Nz = p.Ns
    dψV = i.fsallast
    dψ = @view dψV[1:δNp]
    V = @view dψV[δNp .+ (1:Nz+1)]
    Vmax = maximum(abs.(extrema(V)))
    δ = @view ψδ[δNp .+ (1:Nz+1)]
    ψ = @view ψδ[1:δNp]
    τf = p.τf
    P = p.P
    q = p.q
    δlf = p.δlf
 
    θ = (p.RSD_RS * exp.((ψ .- p.RSf0) ./ p.RSb)) / p.RSV0  # Invert ψ for θ.

    # data for global.dat file
    uu = Array{Float64}(undef, 1, 3)
    uu[1] = t
    uu[2] = NaNMath.log10(Vmax)
    uu[3] = moment_density_rate(V, μshear, dz)
    open(path_to_global, "a") do io
      writedlm(io, uu)
    end
    
    if mod(ctr[], p.save_stride_fields) == 0 || t == (sim_years ./ 31556926)
      vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = NaNMath.log10(Vmax)
      vv[3:end] = δ[flt_loc_indices]
      open(path_to_slip, "a") do io
        writedlm(io, vv)
      end

      stations = ["-15", "+00", "+05", "+10", "+15", "+25", "+35", "+50", "+75"]
      
      for i = 1:length(station_indices)
        ww = Array{Float64}(undef, 1, 7)
        ww[1] = t
        ww[2] = δ[station_indices[i]]
        ww[3] = NaNMath.(V[station_indices[i]])
        ww[4] = τf[station_indices[i]]
        ww[5] = P[station_indices[i]]
        ww[6] = q[station_indices[i]]
        ww[7] = NaNMath.log10(θ[station_indices[i]-δlf+1])

        XXX = pth * "fltst_strk"*stations[i]*".txt"
        open(XXX, "a") do io
            writedlm(io, ww)
        end
      end
    end
  
      
  
    global ctr[] += 1
    @show ctr[]
  

  end
     Vmax

  
  
end


    
function create_text_files_BP6(pth, flt_loc, flt_loc_indices, stations, station_indices, t, RSVinit, δ, τz0, θ, δlf, P, q, μshear)

  path_to_global = pth * "global.dat"
  path_to_slip = pth * "slip.dat"
  # global.dat includes time series of maximum amplitude of slip rates, and moment density rates
  uu = Array{Float64}(undef, 1, 3)
  uu[1] = t
  uu[2] = NaNMath.log10(RSVinit)   # V = V_init everywhere
  uu[3] = μshear * RSVinit * 40 * 1e12 # constants come out of integral, int(dz) = length of RS domain = 40 km
  open(path_to_global, "w") do io
    # write(io, "# problem=SEAS Benchmark BP6-A\n")  # aging law
    write(io, "# problem=SEAS Benchmark BP6-S\n")  # slip law
    write(io, "# code=Thrase\n")
    write(io, "# modeler=J. Marcum\n")
    write(io, "# date=2022/10/20\n")
    write(io, "# element size=100 m\n")
    write(io, "# location=frictional domain\n")
    write(io, "# Column #1 = Time (s)\n")
    write(io, "# Column #2 = Max slip rate (log10 m/s)\n")
    write(io, "# Column #3 = Moment density rate (N/s)\n")
    write(io, "t max_slip_rate moment_rate\n")
    writedlm(io, uu)
  end
  
  # slip.dat is a file that stores time, max(V) and slip at all the stations:
  open(path_to_slip, "w") do io
    write(io,"0.0 0.0 ")
      for i in 1:length(flt_loc)
        write(io,"$(flt_loc[i]) ")
      end
        write(io,"\n")
    end
  
  #write out initial data into devol.txt:
  vv = Array{Float64}(undef, 1, 2+length(flt_loc))
    vv[1] = t
    vv[2] = NaNMath.log10(RSVinit)
    vv[3:end] = δ[flt_loc_indices]
    open(path_to_slip, "a") do io
        writedlm(io, vv)
    end

  # write out initial data into station files:

  # fltst_dpXXX.txt is a file that stores time and time-series of slip, log10(slip_rate), 
  # shear_stress and log10(state) at depth of z = XXX km, where XXX is each of the fault station depths.
  # First we write out initial data into each fltst_dpXXX.txt:

  stations = ["-15", "+00", "+05", "+10", "+15", "+25", "+35", "+50", "+75"]
  for n = 1:length(stations)
    XXX = pth * "fltst_strk"*stations[n]*".txt"
    ww = Array{Float64}(undef, 1, 7)
    ww[1] = t
    ww[2] = δ[station_indices[n]]
    ww[3] = NaNMath.log10(RSVinit)
    ww[4] = τz0
    ww[5] = P[station_indices[n]]
    ww[6] = q[station_indices[n]]
    ww[7] = NaNMath.log10(θ[station_indices[n]-δlf+1])  # subtract off number of points outside RS region?
    open(XXX, "w") do io
      # write(io, "# problem=SEAS Benchmark BP6-A\n")  # aging law
      write(io, "# problem=SEAS Benchmark BP6-S\n")  # slip law
      write(io, "# code=Thrase\n")
      write(io, "# modeler=J. Marcum\n")
      write(io, "# date=2022/10/16\n")
      write(io, "# element size=100 m\n")
      write(io, "# location=on fault, z = "*string(parse(Int64, stations[n])/10)*" km\n")
      write(io, "# Lz = 40 km\n")
      write(io, "t slip slip_rate shear_stress pore_pressure darcy_vel state\n")

      writedlm(io, ww)
    end
  end

end

