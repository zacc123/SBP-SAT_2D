# SBP-SAT_2D

## High Level Goal:
Repo to keep track of 2D SBP-SAT Work for eventual use in the Thrase.jl project from Prof. Brittany Erickson.

## Project Description:
Implement the SBP-SAT methods for the 2D elasticity equations in (Erickson + Dunham, 2014), (Erickson, et al. 2020), (Almquist + Dunham 2020), (Kozdon, et al. 2021) among others (See links to papers below : ) ) using Variable Coefficients, coordinate transforms, and all that fun stuff. Ultimate goal is to take what's implemented here and use it in the [Thrase.jl](https://github.com/Thrase/Thrase.jl) code base for SEAS work.

## Current Progress:
As of 7-10-2025,
- Some peformance changes made to 2D Variable Coefficient code now "2D_Wave_3.1_VarCoeffs_Faster.jl"
    - **This is now the one to run for any experiments**
- Shifted AFC operators and Normal FC operators to come from get_ops_draft_1.1.jl (need to rename tho haha)
- Changes tested with the following settings:
    - p=2,4,6 for $y(r, s) = s, x(r, s) = 2r$ on grid: $y \in (-1, 1), x \in (-2, 2)$
    - Converged to 2, 3, and 4 ~10 min, 10min, and 45 min respectively for (81x81) meshes, on $0 \leq t \leq 2, dt = 1e-4$

    
As of 7-7-2025, 
- MMS updated so $u(x, y, t) = sin(C(x+y) - t)$
- AFC operators implemented for p=2,4,6, ... whatever is implemented in the diag_sbp.jl
- Some Speed Updates (Tho minor : / ) with the most current 2D_Wave_3.1_VarCoeffs_Faster.jl
    - Made some minor memory management to reduce total num of allocations per @time macro by 10%
    - Need to look into better ways to identify bottlenecks and see what I can fix
## TO-DO:
- [x] implement the AFC and standard operators (Mattson 2012), (Erickson et al 2020) for p > 2
- [x] implement MMS for u(x, y, t) = sin(c(x + y) - t) to better match a wave solution
- [ ] Automated Unit Testing for metric, operator, and other identities in the papers

## Papers:
- Erickson, B. A., Kozdon, J. E., and Harvey, T. (2022), A non-stiff summation-by-parts finite difference method for the wave equation in second order form: Characteristic boundary conditions and nonlinear interfaces, Journal of Scientific Computing, doi: 10.1007/s10915-022-01961-1.
- Kozdon, J. E., Erickson, B. A., and Wilcox, L. C. (2020), Hybridized summation-by-parts finite difference methods, Journal of Scientific Computing, doi: 10.1007/s10915-021-01448-5.
- Erickson, B. A. and Dunham, E. M. (2014), An efficient numerical method for earthquake cycles in heterogeneous media: Alternating sub-basin and surface-rupturing events on faults crossing a sedimentary basin, Journal of Geophysical Research, doi:10.1002/2013JB010614.
