# `m_diffusion.f90`

## Overview

The `m_diffusion` module provides high-level solver routines for time-dependent diffusion problems, typically of the form $\frac{\partial \phi}{\partial t} = \nabla \cdot (D \nabla \phi) + S$, where $D$ is the diffusion coefficient (tensor) and $S$ is a source term (though source terms are implicitly handled by the RHS construction). This module transforms the time-discretized diffusion equation into a Helmholtz-type equation, which is then solved using the multigrid capabilities of the `m_octree_mg` library. It supports constant, variable isotropic, and variable anisotropic diffusion coefficients by leveraging specialized Helmholtz solver modules (`m_helmholtz`, `m_vhelmholtz`, `m_ahelmholtz`). Both first-order (e.g., backward Euler) and second-order (e.g., Crank-Nicolson) time discretizations are supported for the diffusion term.

## Key Components

### Modules

- **`m_diffusion`:** Contains public subroutines for solving different types of diffusion equations.

### Functions/Subroutines

The module provides three main public solver routines, each tailored to a different type of diffusion coefficient:

- **`diffusion_solve(mg, dt, diffusion_coeff, order, max_res)`:**
  - **Description:** Solves a diffusion equation where the diffusion coefficient `diffusion_coeff` is constant throughout the domain. It internally sets `mg%operator_type = mg_helmholtz` and uses the `m_helmholtz` module to perform the solution.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure. The solution $\phi^t$ is expected in `mg%boxes(:)%cc(..., mg_iphi)` and is overwritten by $\phi^{t+dt}$.
    - `dt (real(dp), intent(in))`: The time step size.
    - `diffusion_coeff (real(dp), intent(in))`: The constant diffusion coefficient.
    - `order (integer, intent(in))`: Time discretization order (1 or 2).
    - `max_res (real(dp), intent(in))`: The desired maximum residual for the solver.

- **`diffusion_solve_vcoeff(mg, dt, order, max_res)`:**
  - **Description:** Solves a diffusion equation with a variable isotropic diffusion coefficient. The coefficient values are expected to be pre-loaded into the `mg_iveps` component of the `mg%boxes(:)%cc` arrays at all relevant grid levels. It internally sets `mg%operator_type = mg_vhelmholtz` and uses the `m_vhelmholtz` module.
  - **Arguments:** (Same as `diffusion_solve`, but `diffusion_coeff` is implicit via `mg_iveps`)
    - `mg (type(mg_t), intent(inout))`
    - `dt (real(dp), intent(in))`
    - `order (integer, intent(in))`
    - `max_res (real(dp), intent(in))`

- **`diffusion_solve_acoeff(mg, dt, order, max_res)`:**
  - **Description:** Solves a diffusion equation with an anisotropic diffusion coefficient. The components of the diffusion tensor are expected to be pre-loaded into `mg_iveps1` (for $D_{xx}$), `mg_iveps2` (for $D_{yy}$), and `mg_iveps3` (for $D_{zz}$) components of `mg%boxes(:)%cc` at all relevant grid levels. It internally sets `mg%operator_type = mg_ahelmholtz` and uses the `m_ahelmholtz` module.
  - **Arguments:** (Same as `diffusion_solve_vcoeff`)
    - `mg (type(mg_t), intent(inout))`
    - `dt (real(dp), intent(in))`
    - `order (integer, intent(in))`
    - `max_res (real(dp), intent(in))`

**Common Solver Logic:**
1.  Set the appropriate `mg%operator_type` (e.g., `mg_helmholtz`).
2.  Call `mg_set_methods(mg)` (from `m_multigrid`) to associate the correct operator and smoother routines from the chosen Helmholtz module with the function pointers in `mg%box_op` and `mg%box_smoother`.
3.  Determine the `lambda` parameter for the Helmholtz equation based on `dt`, `diffusion_coeff` (if applicable), and `order`.
4.  Construct the right-hand side (RHS) of the Helmholtz equation using the private `set_rhs` subroutine. This involves terms from $\phi^t$. For `order=2`, it also includes a term $L_{op}(\phi^t)$, where $L_{op}$ is the spatial diffusion operator part, effectively moving known terms to the RHS.
5.  Execute a Full Multigrid (FMG) cycle using `mg_fas_fmg`.
6.  Iteratively call V-cycles (`mg_fas_vcycle`) until the residual `res` is less than or equal to `max_res`, or a maximum iteration count (`max_its = 10`) is reached.
7.  If convergence is not achieved, an error is raised.

### Private Helper Subroutines

- **`set_rhs(mg, f1, f2)`:**
  - **Description:** This subroutine computes the right-hand side for the Helmholtz-like equation that results from the time discretization of the diffusion equation.
    For a scheme like $(\alpha I - \beta D_0 \Delta t \nabla^2) \phi^{n+1} = \gamma \phi^n + \delta D_0 \Delta t \nabla^2 \phi^n + S'$, this routine effectively calculates $\text{RHS} = \text{f1} \cdot \phi^n + \text{f2} \cdot (\text{previous content of } mg\_irhs)$.
    The `mg_irhs` field in `mg%boxes(:)%cc` is overwritten with this new RHS.
    - If `order = 1`: $\text{RHS} = (1/(\text{dt} \cdot D)) \phi^n$ (for constant D, simplified). `f1 = 1/(\text{dt} \cdot D)`, `f2 = 0`.
    - If `order = 2`: $\text{RHS} = (2/(\text{dt} \cdot D)) \phi^n + \nabla^2 \phi^n$. `f1 = 2/(\text{dt} \cdot D)`, `f2 = -1` (since $\nabla^2 \phi^n$ was stored in `mg_irhs` by `mg_apply_op` with a negative sign from the Helmholtz convention $L\phi = (\nabla^2 - \lambda)\phi$).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `f1 (real(dp), intent(in))`: Factor multiplying the current solution $\phi^t$ (stored in `mg_iphi`).
    - `f2 (real(dp), intent(in))`: Factor multiplying the previous content of `mg_irhs` (which for `order=2` holds the result of applying the spatial operator to $\phi^t$).

## Important Variables/Constants

- **`max_its (integer, parameter, local to each public solver)`:** Maximum number of V-cycles allowed after the initial FMG cycle. Hardcoded to `10`.
- **Variable Interpretation:**
  - `mg%boxes(:)%cc(..., mg_iphi)`: Input is $\phi^t$, output is $\phi^{t+dt}$.
  - `mg%boxes(:)%cc(..., mg_irhs)`: Used to construct the RHS of the Helmholtz solver.
  - `mg%boxes(:)%cc(..., mg_iveps*)`: Stores diffusion coefficients for variable/anisotropic cases.
- **`dt (real(dp))`**: Timestep. Critical for stability and accuracy, and for setting Helmholtz $\lambda$.
- **`order (integer)`**: Time discretization order (1 or 2).

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)` and its components (e.g., `mg_iphi`, `mg_irhs`, `mg_iveps`, `operator_type`), and various constants.
  - **`m_multigrid`:** Crucial for `mg_set_methods` (to link the correct Helmholtz operator/smoother), `mg_apply_op` (used in `order=2` schemes to compute $L\phi^t$), `mg_fas_fmg`, and `mg_fas_vcycle` (the core multigrid solution algorithms).
  - **`m_helmholtz`:** Provides `helmholtz_set_lambda`; its operator/smoother are used by `diffusion_solve`.
  - **`m_vhelmholtz`:** Provides `vhelmholtz_set_lambda`; its operator/smoother are used by `diffusion_solve_vcoeff`.
  - **`m_ahelmholtz`:** Provides `ahelmholtz_set_lambda`; its operator/smoother are used by `diffusion_solve_acoeff`.
- **External Libraries:**
  - **`cpp_macros.h`:** Used for the `DTIMES` macro in `set_rhs` for dimension-agnostic array indexing.
- **Interactions with Other Components:**
  - This module serves as a high-level application driver that orchestrates the use of underlying multigrid solvers for a specific physical problem (diffusion).
  - **Calling Code:** The user must:
    1.  Initialize the `mg_t` structure (grid, MPI setup, boundary conditions via `mg%bc`).
    2.  Store the initial condition $\phi^t$ in `mg%boxes(:)%cc(..., mg_iphi)`.
    3.  For `diffusion_solve_vcoeff` or `diffusion_solve_acoeff`, ensure the diffusion coefficients are correctly populated in the `mg_iveps*` fields of `mg%boxes(:)%cc` *on all grid levels* that will be used by the solver. This often requires defining how coefficients are restricted to coarser grids.
  - **Output:** The new solution $\phi^{t+dt}$ is stored in `mg%boxes(:)%cc(..., mg_iphi)`.

## Usage Examples

```fortran
! Conceptual example: Solving a time step of diffusion
! with a constant diffusion coefficient.

use m_data_structures
use m_diffusion
! Assume 'mg_state' is a type(mg_t) variable that has been fully initialized:
! - MPI setup done (e.g., via mg_comm_init)
! - Tree built (e.g., via mg_build_rectangle)
! - Storage allocated (e.g., via mg_allocate_storage)
! - Boundary conditions for mg_iphi are set in mg_state%bc
! - Initial solution phi_t is loaded into mg_state%boxes(:)%cc(..., mg_iphi)

implicit none

type(mg_t) :: mg_state
real(dp)   :: time_step_val
real(dp)   :: const_diffusion_coeff
integer    :: discretization_order_val
real(dp)   :: convergence_tolerance

! Parameters for the diffusion solve
time_step_val            = 0.001_dp
const_diffusion_coeff    = 0.1_dp
discretization_order_val = 2  ! Use a second-order time scheme (e.g., Crank-Nicolson like)
convergence_tolerance    = 1.0e-8_dp

! Call the diffusion solver for one time step
call diffusion_solve(mg_state, time_step_val, const_diffusion_coeff, &
                     discretization_order_val, convergence_tolerance)

! The updated solution phi_(t+dt) is now in mg_state%boxes(:)%cc(..., mg_iphi)
print *, "Diffusion step completed. Solution updated."

! To solve with a variable coefficient (assuming D is in mg_iveps):
! call diffusion_solve_vcoeff(mg_state, time_step_val, &
!                             discretization_order_val, convergence_tolerance)

! To solve with an anisotropic coefficient (assuming Dxx, Dyy, Dzz are in mg_iveps1/2/3):
! call diffusion_solve_acoeff(mg_state, time_step_val, &
!                             discretization_order_val, convergence_tolerance)
```
