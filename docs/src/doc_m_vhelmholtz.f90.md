# `m_vhelmholtz.f90`

## Overview

The `m_vhelmholtz` module implements procedures for the solution of a variable-coefficient Helmholtz equation using the `m_octree_mg` multigrid library. The general form of the equation solved is $\nabla \cdot (D(\mathbf{x}) \nabla \phi) - \lambda \phi = f$. In this equation:
- $\phi$ is the unknown field to be solved.
- $D(\mathbf{x})$ is a spatially varying isotropic diffusion coefficient. This coefficient must be provided by the user and stored as a cell-centered quantity in the `mg_iveps` variable component of the `mg%boxes(:)%cc` data arrays for each grid box.
- $\lambda$ is a scalar constant (non-negative).
- $f$ is the source term.

This module provides the discretized form of this operator and a corresponding Gauss-Seidel based smoother, which are then utilized by the main multigrid solver routines. The discretization uses harmonic averaging for the diffusion coefficient $D(\mathbf{x})$ at cell faces.

## Key Components

### Modules

- **`m_vhelmholtz`:** Encapsulates the variable-coefficient Helmholtz operator, smoother, and configuration routines.

### Public Module Variable

- **`vhelmholtz_lambda (real(dp), public, protected)`:**
  - **Description:** Stores the global constant value for the $\lambda$ term in the Helmholtz equation. This value must be non-negative ($\ge 0$). It is initialized to `0.0_dp`. If $\lambda = 0$, the equation becomes a variable-coefficient Poisson equation $\nabla \cdot (D(\mathbf{x}) \nabla \phi) = f$.
  - **Accessed by:** `vhelmholtz_set_lambda` (to set) and the private `box_gs_vhelmh` and `box_vhelmh` routines (to use).

### Functions/Subroutines

- **`vhelmholtz_set_methods(mg)`:**
  - **Description:** Configures the provided multigrid structure (`mg` of type `mg_t`) to use the variable-coefficient Helmholtz operator and smoother defined in this module.
    - It assigns the internal `box_vhelmh` subroutine to `mg%box_op` (operator application) and `box_gs_vhelmh` to `mg%box_smoother`.
    - It checks if `mg%n_extra_vars` is at least 1 (required for storing the diffusion coefficient $D(\mathbf{x})$ in `mg_iveps`). If `mg` is already allocated and `n_extra_vars` is 0, it stops with an error. Otherwise, it ensures `n_extra_vars` is at least 1.
    - It sets `mg%subtract_mean` to `.false.`.
    - Crucially, it sets the boundary condition type for the `mg_iveps` variable to `mg_bc_neumann` with a value of `0.0_dp`. This ensures that the diffusion coefficient $D$ can be correctly accessed in ghost cells when applying stencils near box boundaries, typically by copying the edge value into the ghost cell.
    - Currently supports `mg_cartesian` geometry and smoother types `mg_smoother_gs` or `mg_smoother_gsrb`.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.

- **`vhelmholtz_set_lambda(lambda)`:**
  - **Description:** Sets the value of the module-level variable `vhelmholtz_lambda`. It enforces that `lambda` must be non-negative; an error is raised if `lambda < 0`.
  - **Arguments:**
    - `lambda (real(dp), intent(in))`: The value for the constant $\lambda$.

### Private Module Procedures (Core Operations)

- **`box_gs_vhelmh(mg, id, nc, redblack_cntr)`:**
  - **Description:** Performs a Gauss-Seidel relaxation sweep (standard or Red-Black, based on `mg%smoother_type`) for the variable-coefficient Helmholtz equation on the data within a single box `id`. It updates `mg_iphi` using values from neighboring cells, the source term `mg_irhs`, the spatially varying diffusion coefficient $D$ (from `mg_iveps`), and `vhelmholtz_lambda`. The diffusion coefficient $D$ at cell faces is computed using harmonic averaging: $D_{face} = 2 D_0 D_{neighbor} / (D_0 + D_{neighbor})$, where $D_0$ is $D$ at the current cell and $D_{neighbor}$ is $D$ at the neighboring cell.
  - **Arguments:** (Standard signature for `mg_box_gsrb`)

- **`box_vhelmh(mg, id, nc, i_out)`:**
  - **Description:** Applies the discretized variable-coefficient Helmholtz operator ($\nabla \cdot (D \nabla \phi) - \lambda \phi$) to the solution variable `phi` (in `mg_iphi`) within box `id`. The result is stored in `cc(..., i_out)`. Like the smoother, it uses harmonic averaging for the diffusion coefficient $D$ at cell faces.
  - **Arguments:** (Standard signature for `mg_box_op`)

## Important Variables/Constants

- **`vhelmholtz_lambda (real(dp))`**: The constant scalar value $\lambda$.
- **Diffusion Coefficient $D(\mathbf{x})$ Storage:** The spatially varying diffusion coefficient is expected to be stored by the user in the `mg_iveps` component of the `mg%boxes(id)%cc` array for every cell in every box.
- **Harmonic Averaging:** The discretization of $\nabla \cdot (D \nabla \phi)$ uses harmonic averaging for $D$ at cell interfaces. For a face between cell $i$ and cell $i+1$, the effective $D_{face}$ is $2 D_i D_{i+1} / (D_i + D_{i+1})$. This is important for maintaining accuracy when $D$ varies sharply.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, `dp`, `NDIM`, variable indices (`mg_iphi`, `mg_irhs`, `mg_iveps`), geometry types (`mg_cartesian`), smoother types (`mg_smoother_gs`, `mg_smoother_gsrb`), and boundary condition types (`mg_bc_neumann`).
- **External Libraries:**
  - **`cpp_macros.h`:** Used for conditional compilation (`#if NDIM == ...`) to provide dimension-specific implementations of stencil operations.
- **Interactions with Other Components:**
  - **`m_multigrid`:** This is the primary consumer module. `vhelmholtz_set_methods` is called (usually via `mg_set_methods` in `m_multigrid` when `mg%operator_type` is `mg_vhelmholtz`) to register `box_vhelmh` as `mg%box_op` and `box_gs_vhelmh` as `mg%box_smoother`. The generic multigrid cycling routines then invoke these.
  - **`m_ghost_cells`:** Before the multigrid routines call the operator or smoother, ghost cells for `mg_iphi`, `mg_irhs`, and critically for `mg_iveps` (the variable diffusion coefficient) must be filled. `vhelmholtz_set_methods` configures Neumann boundary conditions for `mg_iveps` to facilitate this.
  - **Application Modules (e.g., `m_diffusion`):** The `diffusion_solve_vcoeff` routine in `m_diffusion` uses `m_vhelmholtz` to solve time-dependent diffusion equations where the diffusion coefficient varies spatially. The application code is responsible for populating `mg_iveps` at all grid levels.

## Usage Examples

```fortran
! Conceptual Example: Setting up to solve a variable-coefficient Helmholtz problem.

use m_data_structures
use m_vhelmholtz
use m_multigrid      ! For the generic mg_set_methods

implicit none

type(mg_t) :: my_problem_config
real(dp)   :: const_lambda_term
integer    :: lvl, i, id, nc
integer    :: iveps_idx ! Should be mg_iveps

! 1. Initialize 'my_problem_config' (grid, MPI, etc.)
!    call comprehensive_mg_initialization(my_problem_config) ! Placeholder
my_problem_config%geometry_type = mg_cartesian
my_problem_config%smoother_type = mg_smoother_gsrb

! 2. Ensure space for the variable coefficient D(x) is requested BEFORE allocation.
!    mg_iveps is typically 5.
iveps_idx = mg_iveps
my_problem_config%n_extra_vars = max(my_problem_config%n_extra_vars, iveps_idx - mg_num_vars + 1)
!    ... (call mg_allocate_storage(my_problem_config) ) ...

! 3. Populate the diffusion coefficient D(x) and the source term f.
!    This must be done for all relevant boxes and levels.
do lvl = my_problem_config%lowest_lvl, my_problem_config%highest_lvl
  nc = my_problem_config%box_size_lvl(lvl)
  do i = 1, size(my_problem_config%lvls(lvl)%my_ids)
    id = my_problem_config%lvls(lvl)%my_ids(i)
    ! Example: Set D(x) = 1.0 + coordinate_x, and f = some_function
    ! my_problem_config%boxes(id)%cc(..., iveps_idx) = 1.0_dp + my_problem_config%boxes(id)%r_min(1) + ...
    ! my_problem_config%boxes(id)%cc(..., mg_irhs) = ...
  end do
end do

! 4. Specify that the problem is a variable-coefficient Helmholtz type.
my_problem_config%operator_type = mg_vhelmholtz

! 5. Set the constant lambda for the Helmholtz term.
const_lambda_term = 1.5_dp
call vhelmholtz_set_lambda(const_lambda_term)

! 6. Associate the vhelmholtz operator and smoother with the mg_t object.
!    This will also set up Neumann BCs for mg_iveps.
call mg_set_methods(my_problem_config)

! Now, 'my_problem_config' is configured. When multigrid routines
! (e.g., mg_fas_fmg) are called, they will use 'box_vhelmh' and 'box_gs_vhelmh'.
! Crucially, ghost cells for mg_iphi, mg_irhs, AND mg_iveps must be filled
! during the multigrid cycles (call mg_fill_ghost_cells_lvl(..., mg_iveps)).
! The Neumann BC setup in vhelmholtz_set_methods helps with mg_iveps ghost cells.
```
