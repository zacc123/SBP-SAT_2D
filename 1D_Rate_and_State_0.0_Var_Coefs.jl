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
include("./RS_function.jl")

# include("./unit_tests.jl")

global const Lx::Int64 = 80
global const AFC::Bool = true
global const C::Float64 = 1
global const year_seconds::Float64 = 60 * 60 * 24 * 365.25
global const ctr = Ref{Int64}(1) 

# RHS above mimics ode_fun fronm thrase, mostly copied and adjusted at this point
# Might have some issues with immutability and params, may need to adjust scalars to vectors

function main()

    (pth, stride_space, stride_time, xc, zc, NX, Nz,
    sim_years, Vp, ρ, cs, σn, RSamin, RSamax, RSb, RSDc,
    RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])
    #########################
    # INITIAL PROBLEM SETUP #
    #########################

    μ = cs^2 / ρ
    η = μ / (2 * cs)

    # Rewrite things how I use it below
    dx = (xc[2] - xc[1]) / NX
    X0, XN = xc

    T0, TN, dT = (0, sim_years * 3600*24*365, 1) # Time grid setup

    X_GRID = X0:dx:XN # store range objects
    T_GRID = T0:dT:TN 

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

    D2, B, T, e = sbp_operators_1D(2, X0, XN, NX, dx, μ, metrics)

    # Get logical coordinates
    dr = 2.0 / NXp
    rc = (-1, 1, dr)
    R_GRID = -1:dr:1
    
    b = zeros(NXp)

    t = 0.0

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

    plot_fault_time_series("slip", pth*"fltst_strk000.txt")
    png("./res/slip.png")
    plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
    png("./res/slip_rate.png")
    
    return nothing
end

main()
# examples of how ot plot times series of shear stress:
    