# SBP-SAT_2D

## High Level Goal:
Repo to keep track of 2D SBP-SAT Work for eventual use in the Thrase.jl project from Prof. Brittany Erickson.

## Project Description:
Implement the SBP-SAT methods for the 2D elasticity equations in (Erickson + Dunham, 2014), (Erickson, et al. 2020), (Almquist + Dunham 2020), (Kozdon, et al. 2021) among others (See links to papers below : ) ) using Variable Coefficients, coordinate transforms, and all that fun stuff. Ultimate goal is to take what's implemented here and use it in the [Thrase.jl](https://github.com/Thrase/Thrase.jl) code base for SEAS work.

## Current Progress:
- As of 7-2-2025, Operators are implemented for Adapted Fully Compatiable case (Alquist and Dunham 2020) for order 2. Current code uses the Method of Manufactured Solution for
    $u(x, y, t) = x^3 + y^3 + t^2$,
    $x \in (-4, -4), y \in (-2, 2)$ with coordinate transformation $x(r, s) = 4r$, $y(r, s) = 2s$
- Code is showing same convergence at p=2, but has not been tested for p=4, 6.

## TO-DO:
- [ ] implement the AFC and standard operators (Mattson 2012), (Erickson et al 2020) for p > 2
- [ ] implement MMS for u(x, y, t) = sin(c(x + y) - t) to better match a wave solution
- [ ] Automated Unit Testing for metric, operator, and other identities in the papers

## Papers:
- Erickson, B. A., Kozdon, J. E., and Harvey, T. (2022), A non-stiff summation-by-parts finite difference method for the wave equation in second order form: Characteristic boundary conditions and nonlinear interfaces, Journal of Scientific Computing, doi: 10.1007/s10915-022-01961-1.
- Kozdon, J. E., Erickson, B. A., and Wilcox, L. C. (2020), Hybridized summation-by-parts finite difference methods, Journal of Scientific Computing, doi: 10.1007/s10915-021-01448-5.
- Erickson, B. A. and Dunham, E. M. (2014), An efficient numerical method for earthquake cycles in heterogeneous media: Alternating sub-basin and surface-rupturing events on faults crossing a sedimentary basin, Journal of Geophysical Research, doi:10.1002/2013JB010614.


