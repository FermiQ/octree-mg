# `m_helmholtz.f90`

## Overview

The `m_helmholtz` module is responsible for implementing the constant-coefficient Helmholtz operator and its associated smoother within the `m_octree_mg` multigrid library. The Helmholtz equation typically takes the form $\nabla^2\phi - \lambda\phi = f$ (or, with a sign convention difference, $\Delta\phi + \lambda\phi = f$), where $\nabla^2$ (or $\Delta$) is the Laplacian operator, $\phi$ is the unknown field, $f$ is the source term, and $\lambda$ is a scalar constant. This module provides the discretized form of this operator and a Gauss-Seidel based smoother tailored for it, which can then be used by the main multigrid solver routines. If $\lambda = 0$, this module effectively solves the Poisson equation.

## Key Components

### Modules

- **`m_helmholtz`:** Encapsulates the Helmholtz operator, smoother, and configuration routines.

### Public Module Variable

- **`helmholtz_lambda (real(dp), public, protected)`:**
  - **Description:** Stores the global constant value for $\lambda$ used in the Helmholtz equation $\nabla^2\phi - \lambda\phi = f$. This value must be non-negative ($\ge 0$). It is initialized to `0.0_dp`.
  - **Accessed by:** `helmholtz_set_lambda` (to set) and the private `box_gs_helmh` and `box_helmh` routines (to use).

### Functions/Subroutines

- **`helmholtz_set_methods(mg)`:**
  - **Description:** Configures the provided multigrid structure (`mg` of type `mg_t`) to use the Helmholtz operator and smoother defined in this module. It achieves this by assigning the internal `box_helmh` subroutine to the `mg%box_op` procedure pointer and the `box_gs_helmh` subroutine to the `mg%box_smoother` procedure pointer. This routine currently supports `mg_cartesian` geometry and smoother types `mg_smoother_gs` (Gauss-Seidel) or `mg_smoother_gsrb` (Gauss-Seidel Red-Black).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.

- **`helmholtz_set_lambda(lambda)`:**
  - **Description:** Sets the value of the module-level variable `helmholtz_lambda`. It checks if the input `lambda` is non-negative; if `lambda < 0`, it stops with an error.
  - **Arguments:**
    - `lambda (real(dp), intent(in))`: The value to be used for $\lambda$ in the Helmholtz equation.

### Private Module Procedures (Core Operations)

- **`box_gs_helmh(mg, id, nc, redblack_cntr)`:**
  - **Description:** This subroutine performs a Gauss-Seidel relaxation sweep (either standard or Red-Black, based on `mg%smoother_type`) for the Helmholtz equation on the data within a single computational box specified by `id`. It updates the solution variable (`mg_iphi`) using values from neighboring cells and the source term (`mg_irhs`). This is the core smoother routine used in multigrid cycles.
  - **Arguments:** (Standard signature for `mg_box_gsrb` procedure pointer)
    - `mg (type(mg_t), intent(inout))`: The multigrid structure.
    - `id (integer, intent(in))`: ID of the box to process.
    - `nc (integer, intent(in))`: Size of the box (number of cells in each direction).
    - `redblack_cntr (integer, intent(in))`: Iteration counter, used for Red-Black ordering to determine which set of cells (red or black) to update.

- **`box_helmh(mg, id, nc, i_out)`:**
  - **Description:** This subroutine applies the discretized constant-coefficient Helmholtz operator ($\nabla^2\phi - \lambda\phi$) to the solution variable `phi` (stored at index `mg_iphi` in the box's `cc` array). The result of the operation is stored in the `i_out` component of the box's `cc` array. This is the core operator routine used, for example, in residual calculations.
  - **Arguments:** (Standard signature for `mg_box_op` procedure pointer)
    - `mg (type(mg_t), intent(inout))`: The multigrid structure.
    - `id (integer, intent(in))`: ID of the box to process.
    - `nc (integer, intent(in))`: Size of the box.
    - `i_out (integer, intent(in))`: Index in the `cc` array where the operator result $L(\phi)$ will be stored.

## Important Variables/Constants

- **`helmholtz_lambda (real(dp))`**: The primary parameter of this module. If `helmholtz_lambda = 0.0_dp`, the operator becomes the Poisson operator.
- **Discretization:** The Laplacian part ($\nabla^2\phi$) is discretized using a standard 7-point stencil in 3D (5-point in 2D, 3-point in 1D) based on second-order central differences.
- **Smoother Factor:** In `box_gs_helmh`, `fac = 1.0_dp / (2 * sum(idr2) + helmholtz_lambda)` represents the inverse of the diagonal element of the discretized Helmholtz operator matrix, used in the Gauss-Seidel update formula. `idr2` contains $1/dr^2$ for each dimension.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, the kind parameter `dp`, `NDIM`, variable indices (`mg_iphi`, `mg_irhs`), geometry types (`mg_cartesian`), and smoother types (`mg_smoother_gs`, `mg_smoother_gsrb`).
- **External Libraries:**
  - **`cpp_macros.h`:** Used for conditional compilation (`#if NDIM == ...`) to provide dimension-specific implementations of the stencil operations in `box_gs_helmh` and `box_helmh`.
- **Interactions with Other Components:**
  - **`m_multigrid`:** This is the primary consumer module. `helmholtz_set_methods` is called (usually via `mg_set_methods` in `m_multigrid` when `mg%operator_type` is `mg_helmholtz`) to register `box_helmh` as `mg%box_op` and `box_gs_helmh` as `mg%box_smoother`. The generic multigrid cycling routines (`mg_fas_fmg`, `mg_fas_vcycle`) then invoke these procedure pointers.
  - **`m_ghost_cells`:** Before the multigrid cycle calls the operator (`box_helmh`) or the smoother (`box_gs_helmh`), routines from `m_ghost_cells` must be called to ensure that ghost cells for the relevant variables (typically `mg_iphi` and `mg_irhs`) are correctly filled with data from neighboring boxes or physical boundary conditions.
  - **Application Modules (e.g., `m_diffusion`):** Higher-level modules that solve specific physical problems (like time-dependent diffusion) might use `m_helmholtz` as a sub-step. For instance, an implicit time step for a diffusion equation $u_t = D \nabla^2 u$ can be rearranged into a Helmholtz equation for $u^{n+1}$. Such modules would call `helmholtz_set_lambda` and `helmholtz_set_methods` (or `mg_set_methods`) appropriately.

## Usage Examples

```fortran
! Conceptual Example: Configuring the multigrid solver for a Helmholtz problem.

use m_data_structures
use m_helmholtz
use m_multigrid      ! For the generic mg_set_methods, though direct call is also possible.

implicit none

type(mg_t) :: my_mg_problem
real(dp)   :: desired_lambda

! 1. Initialize the multigrid structure (details omitted for brevity)
!    This includes defining the grid, MPI settings, etc.
!    call setup_basic_mg_structure(my_mg_problem)
my_mg_problem%geometry_type = mg_cartesian
my_mg_problem%smoother_type = mg_smoother_gsrb ! e.g., Red-Black Gauss-Seidel

! 2. Specify that the problem is a Helmholtz type
my_mg_problem%operator_type = mg_helmholtz

! 3. Set the Helmholtz lambda parameter
desired_lambda = 5.0_dp
call helmholtz_set_lambda(desired_lambda)

! 4. Associate the Helmholtz operator and smoother with the mg_t object.
!    This can be done by directly calling the module's set_methods:
!    call helmholtz_set_methods(my_mg_problem)
!    Or, more commonly, by calling the generic dispatcher in m_multigrid
!    if mg_operator_type is already set to mg_helmholtz:
call mg_set_methods(my_mg_problem)

! Now, 'my_mg_problem' is configured. When multigrid routines like
! mg_fas_fmg(my_mg_problem, ...) or mg_smooth(my_mg_problem, ...) are called,
! they will use the 'box_helmh' for operator applications and
! 'box_gs_helmh' for smoothing, with the specified 'helmholtz_lambda'.

! Inside a multigrid cycle (conceptual flow):
! integer :: some_box_id, cells_in_box, target_var_index, rb_counter
! ...
! ! Ghost cells for mg_iphi and mg_irhs would be filled here for relevant levels
! ! call mg_fill_ghost_cells_lvl(my_mg_problem, current_level, mg_iphi)
! ! call mg_fill_ghost_cells_lvl(my_mg_problem, current_level, mg_irhs)
!
! ! Operator application (e.g. to compute residual L(phi) - f )
! ! This call would resolve to box_helmh(...)
! call my_mg_problem%box_op(my_mg_problem, some_box_id, cells_in_box, target_var_index)
!
! ! Smoothing step
! ! This call would resolve to box_gs_helmh(...)
! call my_mg_problem%box_smoother(my_mg_problem, some_box_id, cells_in_box, rb_counter)
! ...
```
