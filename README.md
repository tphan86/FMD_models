# Deterministic and stochastic foot-and-mouth (FMD) models

The code was run using the version of MATLAB R2023.

## File descriptions

* `model_setup.mlx`: description of deterministic and stochastic FMD models with tables of variables and parameters
* `g.m`: local function for convex incidence function
* `f_drift.m`: local function for the right-hand side of the deterministic model
* `g1_diff.m`: local function for the diffusion part of exposed population
* `g2_diff.m`: local function for the diffusion part of infectious population
* `RK_stochastic_FMD.m`: local function of stochastic Runge-Kutta algorithm (strong order 1)
* `lyapunov_exponents.m`: local function for computing Lyapunov Exponents of each stochastic solution component
* `threshold.m`: local function for computing the threshold to determine extinction and persistence of the stochastic model
* `jacobian_matrix.m`: local function for computing the Jacobian matrix of the deterministic model at an equilibrium
* `equilibria_ode.m`: local function to compute and check stability of all equilibria of the deterministic model
* `check_stability_sde.m`: local function to check stability of equilibrium states of the stochastic model
* `figure1.m`: produce Figure 2 in the manuscript
* `figure2.m`: produce Figure 3 in the manuscript
* `figure3.m`: produce Figure 4 in the manuscript
* `figure4.m`: produce Figure 5 in the manuscript
* `figure5.m`: produce Figure 6 in the manuscript
* `figure6.m`: produce Figure 7 in the manuscript 
