# Now matches thrase minus z comps
struct ODE_params{a, b, c, d}
    reject_step::a
    sim_years
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

        return nothing
end