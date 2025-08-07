#=
Proof of concept for 1D Rate and State friction problem for scalar elasticity equations
    Using 1D Variable Coefficient SBP-SAT operators from (Alquist and Dunham 2020)
    Using problem setup and time stepping method from (Erickson and Dunham 2014)
=#


# Useful imports
using LinearAlgebra
using Plots
using SparseArrays
using BenchmarkTools
using DifferentialEquations
using NaNMath

# Local files for this project
include("./methods.jl") # Pull in Euler and RK2 Explicit Methods
include("./helper_functions.jl") # Some nice stack, unpack, etc methods
include("./sbp_ops.jl") # We'll actually build our ops
include("./plotting_2D.jl") # Moved all plots to here
include("./convergence_testing.jl") # Moved Convergence here

include("./get_ops_draft_1.1.jl") # Add in the get ops now
include("./coordinate_transform.jl") # Add in coordinate tranforms
include("./1D_ops.jl")

# include("./unit_tests.jl")

global const Lx::Int64 = 80
global const AFC::Bool = true
global const C::Float64 = 1
global const year_seconds::Float64 = 60 * 60 * 24 * 365.25
global const ctr = Ref{Int64}(1) 



# RHS above mimics ode_fun fronm thrase, mostly copied and adjusted at this point
# Might have some issues with immutability and params, may need to adjust scalars to vectors
pth = "./res/"
function main()

    #########################
    # INITIAL PROBLEM SETUP #
    #########################

    ### All Param setting from Thrase Stripped driver
    # how many years to simulate
    sim_years = 2000


    pth = "./res/"
    #
    # loading rate
    Vp = 1e-9
    #
    # elastcity parameters
    ρ = 2.670
    cs = 3.464
    σn = 50.0

    μ = cs^2 / ρ
    η = μ / (2 * cs)

    #
    # rate-and-state friction parameters
    RSamin = 0.01
    RSamax = 0.025
    RSb = 0.015
    RSDc = 0.008
    RSf0 = 0.6
    RSV0 = 1e-6
    RSVinit = 1e-9
    RSH1 = 15
    RSH2 = 18
    RSWf = 40

    SBPp   = 2

    stride_space = 1
    # write-out every "stride_time" time steps
    stride_time = 20
    dx = 10
    # start with grid setup
    X0, XN, dX = (0, Lx, dx) # Physical Grid size 
                              # Remember that logical space goes to x \in (-1, 1)

    T0, TN, dT = (0, sim_years * 3600*24*365, 1) # Time grid setup

    X_GRID = X0:dX:XN # store range objects
    T_GRID = T0:dT:TN 

    NX = ((XN - X0) / dX ) #  Keeping number of nodes with paper conventions
    NX = Int(NX)
    NXp = NX + 1


    # Get SBP related things out of the way
    

    # Setup the coordinate transform 
    function x_to_r(r, s)
        # Scale assuming physical is 0 -> Lx
        # shift to -1, 1
        i, j = size(r)
        x = zeros(i, j)

        # scaling factor
        alpha = XN - X0 / 2 # Scale
        beta = (X0 + 1) * alpha # Shift
        
        # set x
        x .= r
        x .*= alpha
        x .+= beta

        # Partial dervs
        xr = zeros(i, j)
        xr .+= alpha

        xs = zeros(size(s))

        return x, xr, xs
    end

    function y_to_s(r, s)
        # Scale assuming physical is 0 -> Lx
        # shift to -1, 1
        i, j = size(s)
        y = s
        yr = zeros(i, j)
        ys = ones(i, j)

        return y, yr, ys
    end
    
    print("\n Create Metrics: ")
    @time metrics = create_metrics_BP6(SBPp, NX, NX, x_to_r, y_to_s) # 

    # confirm that metrics is working as expected
    @assert metrics.rx[1:NXp] == 1 / (XN - X0 / 2) .* ones(NXp) # make sure that partials are working out
    @assert metrics.sy[1:NXp] == ones(NXp)

    # Next section until # -> will be set in 2D then moved to 1D
    # Drr[1:Nrp, 1:Nrp] will be the only used operator for 1D
    # Get SBP Operators for 2nd Derv
    print("\n Get Ops Almquist: ")
    @time JH, D, H = get_operators_BP6(SBPp, NX, NX, μ, Lx, Lx; metrics=metrics, afc=AFC)

    # ->

    # Get 1st Derv Operators for SAT Terms and Traction
    print("\n")
    Dr_tmp, _, _, _ =  diagonal_sbp_D1(SBPp, NX; xc = (-1, 1))

    Dr_sat = spzeros(NXp, NXp) # Set normal for traction
    Dr_sat[1, :] = Dr_tmp[1, :]
    Dr_sat[end, :] = Dr_tmp[end, :]

    D1s = (Dr_tmp, Dr_sat)
    
    @assert issparse(Dr_sat) && issparse(Dr_tmp)

    # Get logical coordinates
    dr = 2.0 / NXp
    rc = (-1, 1, dr)
    R_GRID = -1:dr:1
    
    # Get remaining SBP Operators from normal, non-var coefficient side
    # TO DO: Condense this down into get ops
    print("Timing for SBP OP Creation:\n")
    @time (_, _, Is, Ir, Hs, Hr, HIs, HIr, BSs, BSr, _, Ef, Er, Es, Ed) = sbp_operators(SBPp, -1, 1, -1, 1, NX, NX, dr, dr)
    
    # Grab the jacobian for only first R Row of S, making it a 1D problem
    J_tmp = zeros(NXp*NXp)
    J = zeros(NXp)
    stack!(J_tmp, metrics.J) # Stack J into same format as R, S
    J .= J_tmp[1:NXp]

    # Get Coef Matrices for SAT Terms
    crr_tmp = spzeros(NXp * NXp, NXp * NXp)
    diagonify(crr_tmp, metrics.crr)
    crr = crr_tmp[1:NXp, 1:NXp]
    @assert issparse(crr)

    D2 = D[1][1:NXp, 1:NXp]
    D2 .*= cs^2
    D2 ./= J
    # D2 = lu(D2)
    @assert issparse(D2)

    # Build the SAT Terms
    Es = Es[1:NXp, 1:NXp] # adjust for 1D
    Ed = Ed[1:NXp, 1:NXp]
    
    # Build Coefficients ahead of time since these will be multiplied by U - Boundary
    alpha_r = -13 / dr
    sat_coef1a = HIr * transpose((alpha_r .* crr) + (crr *  Dr_sat)) * Es ./ J
    sat_coef1b = HIr * transpose((alpha_r .* crr) + (crr *  Dr_sat)) * Ed ./ J

    #print(Array(sat_coef1a))
    # Account for QD case where 0 = D2u + SATu + f
    SAT = sat_coef1a + sat_coef1b
    D2 += SAT

    b = zeros(NXp)
    B = [sat_coef1a, sat_coef1b]
    t = 0.0
    T = [μ .* Dr_tmp] # traction term with mu added!
    e = (Es, Ed)

    δ = zeros(NXp)

    bdry_vec_strip!(b, B, X_GRID, δ ./ 2, (t .* Vp./2)*ones(NXp), Lx)

    u = D2 \ b
    # That should be it as far as SBP SAT Prep

    # initialize change in shear stress due to quasi-static deformation
    Δτ = zeros(1)

    RSa = zeros(1)
    RSa[1] = RSamin - (RSamin - RSamax) * min(1, max(0, (RSH1 - 0)/(RSH1 - RSH2)))
    
    # Set pre-stress according to benchmark description
    τ0 = σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                    exp((RSf0 + RSb * log(RSV0 / RSVinit)) /
                                        RSamax)) + η * RSVinit

    # Set initial state variable according to benchmark
    θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
        sinh.((τ0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)

    # Initialize psi version of state variable
    ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)
        
    # Set initial condition for index 1 DAE - this is a stacked vector of psi, followed by slip
    ψδ = zeros(2)  #because length(ψ) = δNp,  length(δ) = Nz+1
    
    
    ψδ[1] = ψ[1]
    
    ψδ[2] = δ[1]
    
    tspan = (T0, TN)
    odeparam = ODE_params( [false], 
                            sim_years,
                            Vp,
                            D2,
                            u,
                            Δτ,
                            τ0*ones(1),
                            b,
                            μ,
                            RSa,
                            RSb,
                            σn,
                            η,
                            RSV0,
                            τ0,
                            RSDc,
                            RSf0,
                            B,
                            X_GRID,
                            T,
                            e,
                            Lx,
                            stride_time # save every save_stride_fields time steps
                        )

                         # Set fault station locations (depths) specified in benchmark
    stations = [0] # km

    # Function that finds the depth-index corresponding to a station location
    function find_station_index(stations, grid_points)
        numstations = length(stations)
        station_ind = zeros(numstations)
        for i in range(1, stop=numstations)
          station_ind[i] = argmin(abs.(grid_points .- stations[i]))
          station_ind[i]
        end
        return Integer.(station_ind)
      end
      
    flt_loc = [1]  # physical stations (units of km)
    flt_loc_indices = find_station_index(flt_loc, ones(1))
    station_indices = find_station_index(stations, ones(1))
    station_strings = ["000"] # "125" corresponds to 12.5 km down dip; these are necessary for writing to files
  
    # Set call-back function so that fields are written to text file after successful time step only.
    cb_fun = SavingCallback((ψδ, t, i) -> write_to_file(pth, ψδ, t, i, ones(1), flt_loc, flt_loc_indices,station_strings, station_indices, odeparam, "BP1_", 0.1 * year_seconds), SavedValues(Float64, Float64))
  

    #alg = Tsit5() #TSIT5 Doesnt preserve my second order convergence : (
    prob = ODEProblem(RHS, ψδ, tspan, odeparam)
    create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, 0, RSVinit, δ, τ0, θ)
    # Solve DAE using Tsit5(), an adaptive Runge-Kutta method
    sol = solve(prob, Tsit5(); dt=0.2,
             abstol = 1e-5, reltol = 1e-5, save_everystep=true, gamma = 0.2,
             internalnorm=(x, _)->norm(x, Inf), callback=cb_fun)

    plot(sol.t, sol[2, :])
    png("./1D_RS_Results/delta.png")

    plot(log.(sol.t), NaNMath.log.(sol[2, :]))
    png("./1D_RS_Results/delta_log.png")

    plot(sol.t, sol[1, :])
    png("./1D_RS_Results/psi.png")

    plot(log.(sol.t), NaNMath.log.(sol[1, :]))
    png("./1D_RS_Results/psi_log.png")
    #@time sol = solve(prob, alg; abstol=1e-10, reltol=1e-10)

    # examples of how ot plot times series of shear stress:
    plot_fault_time_series("slip", pth*"fltst_strk000.txt")
    png("./res/slip.png")
    plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
    png("./res/slip_rate.png")
    
    return nothing
end

# Now matches thrase minus z comps
struct ODE_params{a, b, c, d}
    reject_step::a
    sim_years::Int64
    Vp::b # array vector
    D2::c # c is sparse array
    u::d
    Δτ
    τf
    b
    μshear
    RSa
    RSb
    σn
    η
    RSV0
    τ0
    RSDc
    RSf0
    B
    x
    T
    e
    Lx    
    save_stride_fields
end

function RHS(dψV, ψδ, params, t) # header now matching Thrase

        # Start with Unpacking
        Vp = params.Vp
        A = params.D2
        u = params.u
        Δτ = params.Δτ
        τf = params.τf
        b = params.b
        μshear = params.μshear
        RSa = params.RSa
        RSb = params.RSb
        σn = params.σn
        η = params.η
        RSV0 = params.RSV0
        τ0 = params.τ0
        RSDc = params.RSDc
        RSf0 = params.RSf0
        B = params.B
        x = params.x 
        T = params.T
        e = params.e
        Lx = params.Lx

        current_time = t ./ 31556926
        print("TIME [YRS] = $(current_time).\n")

        ψ  = ψδ[1]
        δ = zeros(size(x))
        δ[1] = ψδ[2]

        

        remote = zeros(size(x))
        remote[end] = (t .* Vp./2)
      
        
        #print("\nDelta before:", size(δ), "\n")
        bdry_vec_strip!(b, B, x, δ ./ 2, remote, Lx)

        # Step 3... Solve for U in the domain
        u[:] = A \ b


        # set up rates of change for  state and slip
        dψ = dψV[1]
        V  = dψV[2]
        
        dψ = 0 # initialize values to 0
        V  = 0 # initialize values to 0
        
        # Update the fault data
        Δτ .= 0
        
        #print("\nSize of Delta Tau:", size(Δτ), '\n')
        Δτ .= computetraction_stripped(T, u, e)
        τf .= τ0 .+ Δτ

        # Do safe-guarded Newton at every node in rate-and-state friction zone in order to solve for slip rate V.
        ψn = ψ
        an = RSa

        τn = (Δτ .+ τ0)[1]
    
        VR = abs(τn / η)
        VL = -VR
        Vn = V
        obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
        (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                    atolx = 1e-9, rtolx = 1e-9)
        V = Vn # update slip rate
        dψV[2] = Vn
        dψV[1] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0) # update aging law

        if  abs(current_time - 2.000) < 0.02
            plot(x, u)
            png("./1D_RS_Results/u_2.png")
            print("Check 2:\n", u,"\n")
        end


        if  abs(current_time - 20.000) < 0.02
            plot(x, u)
            png("./1D_RS_Results/u_20.png")
            print("Check 20:\n", u,"\n")
        end

        if  abs(current_time - 200.000) < 4 && current_time < 200
            plot(x, u)
            png("./1D_RS_Results/u_196.png")
            print("Check 196:\n", u,"\n")
        end
        if  abs(current_time - 200.000) < 1 && current_time > 200
            plot(x, u)
            print("Check 201:\n", u,"\n")
            png("./1D_RS_Results/u_200.png")
        end

        if  abs(current_time - 2000.000) < 4 && current_time > 200
            plot(x, u)
            print("Check 200:\n", u,"\n")
            png("./1D_RS_Results/u_2000.png")
        end
        return nothing
    end


main()
# examples of how ot plot times series of shear stress:
    plot_fault_time_series("slip", pth*"fltst_strk000.txt")
    png("./res/slip.png")
    plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
    png("./res/slip_rate.png")