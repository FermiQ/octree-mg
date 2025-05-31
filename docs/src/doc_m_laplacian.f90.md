# `m_laplacian.f90`

## Overview

The `m_laplacian` module is dedicated to implementing the Laplacian operator ($\nabla^2 \phi$ or $\Delta \phi$) and associated smoother routines for use within the `m_octree_mg` multigrid library. The Laplacian is a fundamental elliptic differential operator that appears in numerous partial differential equations (PDEs), most prominently in Poisson's equation ($\nabla^2 \phi = f$) and as the spatial component of heat (diffusion) and wave equations. This module provides discretizations for Cartesian coordinates (in 1D, 2D, and 3D) and for 2D Cylindrical coordinates $(r,z)$.

## Key Components

### Modules

- **`m_laplacian`:** Contains the subroutines for setting up and applying the Laplacian operator and its smoother.

### Functions/Subroutines

- **`laplacian_set_methods(mg)`:**
  - **Description:** This public subroutine configures the multigrid structure `mg` (of type `mg_t`) to use the appropriate Laplacian operator and smoother based on the `mg%geometry_type`.
    - If `mg%geometry_type` is `mg_cartesian`:
      - `mg%box_op` (the operator application procedure pointer) is set to `box_lpl`.
      - `mg%box_smoother` (the smoother procedure pointer) is set to `box_gs_lpl` (supporting both standard Gauss-Seidel and Red-Black Gauss-Seidel via `mg%smoother_type`).
    - If `mg%geometry_type` is `mg_cylindrical` (and `NDIM == 2`):
      - `mg%box_op` is set to `box_clpl`.
      - `mg%box_smoother` is set to `box_gs_clpl`.
    - For fully periodic problems in Cartesian coordinates (`all(mg%periodic)` is true), it sets `mg%subtract_mean = .true.`. This flag indicates to the solver that the mean of the solution and RHS should be handled appropriately to ensure a unique solution.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.

### Private Module Procedures (Core Operations)

**For Cartesian Coordinates:**
- **`box_lpl(mg, id, nc, i_out)`:**
  - **Description:** Applies the discretized Cartesian Laplacian operator ($\nabla^2\phi$) to the variable stored at index `mg_iphi` within the computational box `id`. The result is stored in the `cc(..., i_out)` component of the box's data array. It uses a standard 3, 5, or 7-point stencil for 1D, 2D, or 3D respectively, based on second-order central differences.
  - **Arguments:** (Standard signature for `mg_box_op`)
    - `mg (type(mg_t), intent(inout))`, `id (integer, intent(in))`, `nc (integer, intent(in))`, `i_out (integer, intent(in))`.

- **`box_gs_lpl(mg, id, nc, redblack_cntr)`:**
  - **Description:** Performs a Gauss-Seidel relaxation sweep (either standard or Red-Black, depending on `mg%smoother_type`) for the Poisson equation ($\nabla^2\phi = f$) on the data within box `id`. It updates `mg%boxes(id)%cc(..., mg_iphi)` using values from neighboring cells and the source term `mg%boxes(id)%cc(..., mg_irhs)`.
  - **Arguments:** (Standard signature for `mg_box_gsrb`)
    - `mg (type(mg_t), intent(inout))`, `id (integer, intent(in))`, `nc (integer, intent(in))`, `redblack_cntr (integer, intent(in))`.

**For 2D Cylindrical Coordinates $(r,z)$ (if `NDIM == 2`):**
- **`box_clpl(mg, id, nc, i_out)`:**
  - **Description:** Applies the 2D cylindrical Laplacian operator $\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial\phi}{\partial r}\right) + \frac{\partial^2\phi}{\partial z^2}$ to the variable `mg_iphi` in box `id`. The radial coordinate $r$ is assumed to be the first dimension, and $z$ the second.
  - **Arguments:** (Standard signature for `mg_box_op`).

- **`box_gs_clpl(mg, id, nc, redblack_cntr)`:**
  - **Description:** Performs Gauss-Seidel relaxation for the 2D cylindrical Poisson equation on box `id`.
  - **Arguments:** (Standard signature for `mg_box_gsrb`).

## Important Variables/Constants

- **Discretization:** The module uses standard second-order central finite difference approximations for all derivatives.
- **Cartesian Smoother Factor:** In `box_gs_lpl`, `fac = 0.5_dp / sum(idr2)` is the reciprocal of the sum of diagonal coefficients ($2 \sum 1/dr_i^2$), used in the Gauss-Seidel update.
- **Cylindrical Coefficients:** `box_clpl` and `box_gs_clpl` use `r_face` (radial position of cell faces) and `r_inv` (inverse of radial cell center positions) to correctly implement the $1/r \cdot \partial/\partial r (r \cdot \partial\phi/\partial r)$ term.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Crucial for the definitions of `type(mg_t)`, `type(mg_box_t)`, the precision kind `dp`, `NDIM`, variable indices (`mg_iphi`, `mg_irhs`), geometry type constants (`mg_cartesian`, `mg_cylindrical`), and smoother type constants (`mg_smoother_gs`, `mg_smoother_gsrb`).
- **External Libraries:**
  - **`cpp_macros.h`:** Used for conditional compilation based on `NDIM` to provide dimension-specific stencil implementations.
- **Interactions with Other Components:**
  - **`m_multigrid`:** The primary interaction point. `laplacian_set_methods` is typically called by `mg_set_methods` (in `m_multigrid`) when `mg%operator_type` is `mg_laplacian`. This registers the appropriate `box_..._lpl` or `box_..._clpl` routines as the active operator (`mg%box_op`) and smoother (`mg%box_smoother`). The main multigrid algorithms (`mg_fas_fmg`, `mg_fas_vcycle`) then use these procedure pointers.
  - **`m_ghost_cells`:** Before any operator application or smoothing sweep performed by this module's routines, ghost cells for the relevant variables (primarily `mg_iphi` for the operator, and both `mg_iphi` and `mg_irhs` for the smoother) must be filled. This is handled by calls to routines in `m_ghost_cells`.
  - **Application Context:** This module is fundamental when solving Poisson's equation or any PDE involving a Laplacian term. Higher-level modules or the main application driver would configure `mg_t` to use this Laplacian operator.

## Usage Examples

```fortran
! Conceptual Example: Setting up the multigrid solver for a Poisson problem
! using the Cartesian Laplacian.

use m_data_structures
use m_laplacian
use m_multigrid      ! For the generic mg_set_methods

implicit none

type(mg_t) :: my_poisson_problem

! 1. Initialize the multigrid structure (grid parameters, MPI setup, etc.)
!    (Details omitted for brevity - this is a complex step involving other modules)
!    call comprehensive_mg_initialization(my_poisson_problem)

! 2. Specify the problem type and geometry for the solver
my_poisson_problem%operator_type = mg_laplacian  ! Indicate a Laplacian-based problem
my_poisson_problem%geometry_type = mg_cartesian  ! Use Cartesian coordinates
my_poisson_problem%smoother_type = mg_smoother_gsrb ! Select Red-Black Gauss-Seidel

! 3. Associate the Laplacian operator and smoother with the mg_t object.
!    This is typically done via the generic mg_set_methods from m_multigrid,
!    which will dispatch to laplacian_set_methods based on operator_type.
call mg_set_methods(my_poisson_problem)

! Now, 'my_poisson_problem' is configured to use the Cartesian Laplacian.
! When multigrid routines (e.g., mg_fas_fmg, mg_fas_vcycle, mg_smooth) are
! called with 'my_poisson_problem', they will internally use:
! - 'box_lpl' for applying the Laplacian operator (e.g., in residual calculation).
! - 'box_gs_lpl' for smoothing steps.

! If solving a 2D problem in cylindrical coordinates:
#if NDIM == 2
! my_poisson_problem%geometry_type = mg_cylindrical
! call mg_set_methods(my_poisson_problem) ! This would set up box_clpl and box_gs_clpl
#endif

! Inside a multigrid cycle (conceptual flow, managed by m_multigrid routines):
! integer :: some_box_id, cells_per_side, output_idx, rb_iter_count
! ...
! ! Ensure ghost cells are up-to-date for the current level and variable mg_iphi
! ! call mg_fill_ghost_cells_lvl(my_poisson_problem, current_lvl, mg_iphi)
!
! ! Apply Laplacian operator (e.g. L(phi) stored in output_idx)
! ! This call would resolve to box_lpl(...) or box_clpl(...)
! call my_poisson_problem%box_op(my_poisson_problem, some_box_id, cells_per_side, output_idx)
!
! ! Apply smoother
! ! This call would resolve to box_gs_lpl(...) or box_gs_clpl(...)
! call my_poisson_problem%box_smoother(my_poisson_problem, some_box_id, cells_per_side, rb_iter_count)
! ...
```
