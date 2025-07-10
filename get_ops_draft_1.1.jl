#= 
My attempt to adjust get_operators() from ops_stripped_BP6
to return D1, D2, BS, S0, SN, A

Questions:
    - Lines 125-131: Why the redundant looking code?
    - Line 183: Why HS when doing Arr, and Hr for Ass
    - What is sim_years in utils.jl line86 supposed to be?
    - Do I need to follow toby's paper exactly for constructing the D2s? Are my H's already good?
    - Am I just building D1, D2, Hs blockwise like the As?
    - Line 232 is CSR vs CRS in the right order??
=#

# Untouched imports and kronecker def
using SparseArrays
using LinearAlgebra

⊗(A,B) = kron(A, B)

# Create metrics does the 2D coordinate transform
function create_metrics_BP6(pm, Nr, Ns,
  xf=(r,s)->(r, ones(size(r)), zeros(size(r))),
  yf=(r,s)->(s, zeros(size(s)), ones(size(s))))
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  # Derivative operators for the metric terms
  @assert pm <= 8
  pp = pm == 6 ? 8 : pm

  r = range(-1, stop=1, length=Nrp)
  s = range(-1, stop=1, length=Nsp)


  # Create the mesh
  r = ones(1, Nsp) ⊗ r
  s = s' ⊗ ones(Nrp)
  (x, xr, xs) = xf(r, s)
  (y, yr, ys) = yf(r, s)
  print(size(r))
  print(size(s))


  J = xr .* ys - xs .* yr
  @assert minimum(J) > 0

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  # variable coefficient matrix components 

  #TODO -> Not actual mu var right??
  crr = J .* (rx .* rx + ry .* ry)
  crs = J .* (sx .* rx + sy .* ry)
  css = J .* (sx .* sx + sy .* sy)

  #
  # Block surface matrices
  #
  (xf1, yf1) = (view(x, 1, :), view(y, 1, :))
  nx1 = -ys[1, :]
  ny1 =  xs[1, :]
  sJ1 = hypot.(nx1, ny1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1

  (xf2, yf2) = (view(x, Nrp, :), view(y, Nrp, :))
  nx2 =  ys[end, :]
  ny2 = -xs[end, :]
  sJ2 = hypot.(nx2, ny2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2

  (xf3, yf3) = (view(x, :, 1), view(y, :, 1))
  nx3 =  yr[:, 1]
  ny3 = -xr[:, 1]
  sJ3 = hypot.(nx3, ny3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3

  (xf4, yf4) = (view(x, :, Nsp), view(y, :, Nsp))
  nx4 = -yr[:, end]
  ny4 =  xr[:, end]
  sJ4 = hypot.(nx4, ny4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4


  (coord = (x,y),
  facecoord = ((xf1, xf2, xf3, xf4), (yf1, yf2, yf3, yf4)),
  crr = crr, css = css, crs = crs,
  J=J,
  sJ = (sJ1, sJ2, sJ3, sJ4),
  nx = (nx1, nx2, nx3, nx4),
  ny = (ny1, ny2, ny3, ny4),
  rx = rx, ry = ry, sx = sx, sy = sy)
end

function get_operators_BP6(p, Nr, Ns, μ, Lx, Lz; metrics=create_metrics_BP6(p,Nr,Ns))

    # In this project, take r = x, s = z (i.e. problem is set up for no coordinate transformation)
    Nrp = Nr + 1  
    Nsp = Ns + 1
    Np = Nrp * Nsp

    # "coefficient" matrices 
    crr = μ * metrics.crr
    css = μ * metrics.css
    crs = metrics.crs
    csr = crs
    J = metrics.J
   
    # Derivative operators for the rest of the computation
    (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
    Qr = Hr * Dr
    QrT = sparse(transpose(Qr))

    (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))  # initially was ` xc = (0,1)` then changed to `xc = (-1,1)` (neither work)
    Qs = Hs * Ds
    QsT = sparse(transpose(Qs))

 
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)
      
    #{{{ Set up the rr derivative matrix
    ISr0 = Array{Int64,1}(undef,0)
    JSr0 = Array{Int64,1}(undef,0)
    VSr0 = Array{Float64,1}(undef,0)
    ISrN = Array{Int64,1}(undef,0)
    JSrN = Array{Int64,1}(undef,0)
    VSrN = Array{Float64,1}(undef,0)

    

    # Making change to Jeremy's code here to also get the HI, H, D2
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, rand(Nrp))
    IArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
    JArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
    VArr = Array{Float64,1}(undef,Nsp * length(Ae.nzval))
    stArr = 0

    ISr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
    JSr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
    VSr0 = Array{Float64,1}(undef,Nsp * length(S0e.nzval))
    stSr0 = 0

    ISrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
    JSrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
    VSrN = Array{Float64,1}(undef,Nsp * length(SNe.nzval))
    stSrN = 0

    # New addition from Zac
    # Pull out dr0 + drN directly from S's 
    # Store as zeros then update
    dr0 = zeros(Nrp)
    drN = zeros(Nrp)
    
    for j = 1:Nsp
      rng = (j-1) * Nrp .+ (1:Nrp)
      (D2rr, S0e, SNe, HIrr, Hrr, Ae, _) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
      
      # Grab d0
      if j == 1
          dr0 = S0e[1, :]  ./ crr[rng] # Make the adjustment so this is the right scale
      end

      # grab dn
      if j == Nsp
          drN = SNe[end, :]  ./ crr[rng]
      end
      # End New addition from Zac

      (Ie, Je, Ve) = findnz(Ae)
      IArr[stArr .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
      JArr[stArr .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
      VArr[stArr .+ (1:length(Ve))] = Hs[j,j] * Ve
      stArr += length(Ve)

      (Ie, Je, Ve) = findnz(S0e)
      ISr0[stSr0 .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
      JSr0[stSr0 .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
      VSr0[stSr0 .+ (1:length(Ve))] =  Hs[j,j] * Ve
      stSr0 += length(Ve)

      (Ie, Je, Ve) = findnz(SNe)
      ISrN[stSrN .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
      JSrN[stSrN .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
      VSrN[stSrN .+ (1:length(Ve))] =  Hs[j,j] * Ve
      stSrN += length(Ve)

    end

    Ãrr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)
    Sr0 = sparse(ISr0[1:stSr0], JSr0[1:stSr0], VSr0[1:stSr0], Np, Np)
    SrN = sparse(ISrN[1:stSrN], JSrN[1:stSrN], VSrN[1:stSrN], Np, Np)
    Sr0T = sparse(JSr0[1:stSr0], ISr0[1:stSr0], VSr0[1:stSr0], Np, Np)
    SrNT = sparse(JSrN[1:stSrN], ISrN[1:stSrN], VSrN[1:stSrN], Np, Np)

    dr0T = dr0'
    drNT = drN'
  
  
    #{{{ Set up the ss derivative matrix
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, rand(Nsp))
    BSs_mini = SNe - S0e
    IAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
    JAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
    VAss = Array{Float64,1}(undef,Nrp * length(Ae.nzval))
    stAss = 0

    ISs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
    JSs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
    VSs0 = Array{Float64,1}(undef,Nrp * length(S0e.nzval))
    stSs0 = 0

    ISsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
    JSsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
    VSsN = Array{Float64,1}(undef,Nrp * length(SNe.nzval))
    stSsN = 0

    ds0 = zeros(Nsp)
    dsN = zeros(Nsp)
 
  for i = 1:Nrp
      rng = i .+ Nrp * (0:Ns)
      (D2, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, css[rng])
      

      # Grab d0
      if i == 1
          ds0 = S0e[1, :] ./ css[rng]
      end

      # grab dn
      if i == Nsp
          dsN = SNe[end, :] ./ css[rng]
        
      end
      # End New addition from Zac

      (Ie, Je, Ve) = findnz(Ae)
      IAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
      JAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
      VAss[stAss .+ (1:length(Ve))] = Hr[i,i] * Ve
      stAss += length(Ve)

      (Ie, Je, Ve) = findnz(S0e)
      ISs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
      JSs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
      VSs0[stSs0 .+ (1:length(Ve))] = Hr[i,i] * Ve
      stSs0 += length(Ve)

      (Ie, Je, Ve) = findnz(SNe)
      ISsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
      JSsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
      VSsN[stSsN .+ (1:length(Ve))] = Hr[i,i] * Ve
      stSsN += length(Ve)
    end

    Ãss = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)
    Ss0 = sparse(ISs0[1:stSs0], JSs0[1:stSs0], VSs0[1:stSs0], Np, Np)
    SsN = sparse(ISsN[1:stSsN], JSsN[1:stSsN], VSsN[1:stSsN], Np, Np)
    Ss0T = sparse(JSs0[1:stSs0], ISs0[1:stSs0], VSs0[1:stSs0], Np, Np)
    SsNT = sparse(JSsN[1:stSsN], ISsN[1:stSsN], VSsN[1:stSsN], Np, Np)

    ds0T = ds0'
    dsNT = dsN' 
    

    #{{{ Set up the sr and rs derivative matrices
    Ãsr = (QsT ⊗ Ir) * sparse(1:length(crs), 1:length(crs), view(crs, :)) * (Is ⊗ Qr)
    Ãrs = (Is ⊗ QrT) * sparse(1:length(csr), 1:length(csr), view(csr, :)) * (Qs ⊗ Ir)
    #}}}

    Ã = Ãrr + Ãss + Ãrs + Ãsr

    Er0 = sparse([1], [1], [1], Nrp, Nrp)
    ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
    Es0 = sparse([1], [1], [1], Nsp, Nsp)
    EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

    er0 = sparse([1  ], [1], [1], Nrp, 1)
    erN = sparse([Nrp], [1], [1], Nrp, 1)
    es0 = sparse([1  ], [1], [1], Nsp, 1)
    esN = sparse([Nsp], [1], [1], Nsp, 1)

    er0T = sparse([1], [1  ], [1], 1, Nrp)
    erNT = sparse([1], [Nrp], [1], 1, Nrp)
    es0T = sparse([1], [1  ], [1], 1, Nsp)
    esNT = sparse([1], [Nsp], [1], 1, Nsp)

    
    # Store coefficient matrices as matrices
    #
    crs0 = sparse(Diagonal(crs[1:Nrp]))
    crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
    csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
    csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

    # Adding 
    crr0 = sparse(Diagonal(crr[1:Nrp]))
    crrN = sparse(Diagonal(crr[Nrp*Ns .+ (1:Nrp)]))
    css0 = sparse(Diagonal(css[1   .+ Nrp*(0:Ns)]))
    cssN = sparse(Diagonal(css[Nrp .+ Nrp*(0:Ns)]))

    #
    # Surface mass matrices
    #
    H1 = Hs
    H1I = HsI

    H2 = Hs
    H2I = HsI

    H3 = Hr
    H3I = HrI

    H4 = Hr
    H4I = HrI

    #
    # Penalty terms
    #
    if p == 2
        l = 2
        β = 0.363636363
        α = 1 / 2
    elseif p == 4
        l = 4
        β = 0.2505765857
        α = 17 / 48
    elseif p == 6
        l = 7
        β = 0.1878687080
        α = 13649 / 43200
    else
        error("unknown order")
    end

    H̃ = Hs ⊗ Hr
    # Throwing these in here to make it match the Jeremy paper when calc Ds
    HsI = Hs \ Is
    HrI = Hr \ Ir

    H̃I = (HsI ⊗ HrI) # Get this to not remake

    # Now get Ds following Kozdon et al.2021
    # For now specify here:
    adpt_fully_comp = 0

    # For now add in a flag
    # adpt fully comp is order 2q on the interior, q - 1 on the boundary
    if adpt_fully_comp == 0                 
      edge_r = (crrN * ErN) - (crr0 * Er0)
      edge_s = (cssN * EsN) - (css0 * Es0)

      # TO DO Fix this to do the correction here
      Drr =  H̃I * (-Ãrr) + (Is ⊗ (HrI * edge_r * Dr)) 
      Dss =  H̃I * (-Ãss) + ((HsI * edge_s * Ds) ⊗ Ir) 

      Drs =  H̃I * (-Ãrs + ( ((crsN * Qs) ⊗ (erN * erNT)) - ((crs0*Qs) ⊗ (er0 * er0T)) ))
      Dsr =  H̃I * (-Ãsr + ( ((esN * esNT) ⊗ (csrN * Qr)) - ((es0 * es0T) ⊗ (csr0 * Qr)) ))

    else
      # Matches all the papers in formation of operators
      # (-A + BS) term basically VV
      Drr =  H̃I * (-Ãrr + ( ((Hs*crrN) ⊗ (erN * drNT)) - ((Hs*crr0) ⊗ (er0 * dr0T)) ))
      Dss =  H̃I * (-Ãss + ( ((esN * dsNT) ⊗ (Hr*cssN)) - ((es0 * ds0T) ⊗ (Hr*css0)) ))

      Drs =  H̃I * (-Ãrs + ( ((crsN * Qs) ⊗ (erN * erNT)) - ((crs0*Qs) ⊗ (er0 * er0T)) ))
      Dsr =  H̃I * (-Ãsr + ( ((esN * esNT) ⊗ (csrN * Qr)) - ((es0 * es0T) ⊗ (csr0 * Qr)) ))
    end

    # Modify the operator to handle the boundary conditions
  
    JH = sparse(1:Np, 1:Np, view(J, :)) * (Hs ⊗ Hr)
    #return (M̃ , F, τ, H̃, HfI_FT)
    (
    JH = JH,
    D = (Drr, Dss, Drs, Dsr),
    H = (Hr, Hs)
    )
end

function bdry_vec_mod_BP6!(g, F, τ, r, s, bc_Dirichlet, bc_Neumann, metrics)
    Nr = length(r)
    Ns = length(s)

    g[:] .= 0
    
    (xf, yf) = metrics.facecoord
    # τ in code is diagonal matrix of penalty parameters

    # # FACE 1 (Neumann) (left):
    # gN = bc_Neumann(1, 0, s, 0, 1)
    # vf = gN  ./ diag(τ[1])
    # g[:] -= F[1] * vf

    # FACE 1 (Dirichlet) (left)
    lf = 1
    vf = bc_Dirichlet(lf, xf[lf], yf[lf])
    g[:] -= F[1] * vf

    # FACE 2 (Dirichlet) (right):
    lf = 2
    vf = bc_Dirichlet(lf, xf[lf], yf[lf])
    g[:] -= F[2] * vf

    # FACE 3 (Dirichlet) (top):
    lf = 3
    vf = bc_Dirichlet(lf, xf[lf], yf[lf])
    g[:] -= F[3] * vf

    # FACE 4 (Dirichlet) (bottom):
    lf = 4
    vf = bc_Dirichlet(lf, xf[lf], yf[lf])
    g[:] -= F[4] * vf 
end


function source_vec_mod_BP6!(g, JH, source, metrics)
  (x, y) = metrics.coord
  g[:] += JH * source(x[:], y[:])
  
end


function computetraction_stripped_BP6(HfI_FT, sJ, τ, lf, u, δ)
    HfI_FT = HfI_FT[lf]
    τf = τ[lf]
    sJ = sJ[lf]
    return (HfI_FT * u + τf * (δ .- δ / 2)) ./ sJ
  end
  


  function check_physical_domain(r_star, s_star, el_r, el_s, Nr, Ns)
    if r_star > 1
      print("error: increase Nr or dx!\n")
      
    end
    
    if r_star ≈ 1
      print("using constant grid spacing in x-direction\n")
      
    end

    if el_r < 2/Nr
        print("error: increase el_r!\n")
        
    end
    
    if s_star > 1
        print("error: increase Ns or dz!\n")
       
    end
    
    if s_star ≈ 1
      print("using constant grid spacing in z-direction")
     
    end

    if el_s < 2/Ns
        print("error: increase el_s!")
        
    end
    return
  end


  function get_tanh_params(xstart, xfinal, L, xistar, el_c)
    
    num = L - xistar*(xfinal - xstart);
    den = tanh((xistar-1)/el_c) + (xistar - 1)*tanh(-1/el_c);
    
    if xistar == 1
        A = 0;
        B = xfinal - xstart;
        C = 0;
    else
        A = num/den;
        B = A*tanh(-1/el_c) + xfinal - xstart;
        C = xfinal - B;
    end
    return (A, B, C)
  end

  export create_metrics_BP6, get_operators_BP6, bdry_vec_mod_BP6!, computetraction_stripped_BP6, check_physical_domain, get_tanh_params