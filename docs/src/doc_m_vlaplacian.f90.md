# `m_vlaplacian.f90`

## Overview

The `m_vlaplacian` module in the `m_octree_mg` library is designed to handle PDEs involving a variable-coefficient Laplacian operator. Specifically, it addresses equations of the form $\nabla \cdot (D(\mathbf{x}) \nabla \phi) = f$, where $D(\mathbf{x})$ is a spatially varying isotropic coefficient (e.g., diffusion coefficient, thermal conductivity, electrical permittivity) and $f$ is a source term. The coefficient $D(\mathbf{x})$ must be provided by the user and is expected to be stored as a cell-centered quantity in the `mg_iveps` variable component of the `mg%boxes(:)%cc` data arrays for each grid box. This module provides the discretized form of the $\nabla \cdot (D(\mathbf{x}) \nabla \cdot)$ operator and a corresponding Gauss-Seidel based smoother, which are essential for the multigrid solution process. The discretization uses harmonic averaging for the coefficient $D(\mathbf{x})$ at cell faces to ensure accuracy, especially when $D$ varies significantly between cells.

This module is structurally very similar to `m_vhelmholtz`, essentially representing the case where the Helmholtz parameter $\lambda=0$.

## Key Components

### Modules

- **`m_vlaplacian`:** Encapsulates the variable-coefficient Laplacian operator, smoother, and configuration routines.

### Functions/Subroutines

- **`vlaplacian_set_methods(mg)`:**
  - **Description:** Configures the provided multigrid structure (`mg` of type `mg_t`) to use the variable-coefficient Laplacian operator and smoother defined in this module.
    - It assigns the internal `box_vlpl` subroutine to `mg%box_op` (operator application) and `box_gs_vlpl` to `mg%box_smoother`.
    - It verifies that `mg%n_extra_vars` is at least 1 (necessary for storing the coefficient $D(\mathbf{x})$ in `mg_iveps`). If `mg` is already allocated and `n_extra_vars` is 0, it issues an error. Otherwise, it ensures `n_extra_vars` is at least 1.
    - It sets `mg%subtract_mean` to `.false.` (mean subtraction is typically for pure Neumann problems with constant coefficients, which is not the primary target here, though specific boundary conditions might still require it at a higher level).
    - It establishes Neumann zero boundary conditions for the `mg_iveps` variable. This is crucial for ensuring that the coefficient $D$ can be correctly accessed in ghost cells during stencil operations near box boundaries.
    - Currently, it supports `mg_cartesian` geometry and smoother types `mg_smoother_gs` (Gauss-Seidel) or `mg_smoother_gsrb` (Gauss-Seidel Red-Black).
    - A commented-out line `! mg%box_prolong => vlpl_prolong` suggests that a custom prolongation method tailored for variable coefficients (`vlpl_prolong`, also commented out in the file) was considered but is not currently active. The default prolongation method from `m_multigrid` would be used.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.

### Private Module Procedures (Core Operations)

- **`box_gs_vlpl(mg, id, nc, redblack_cntr)`:**
  - **Description:** This subroutine performs a Gauss-Seidel relaxation sweep (either standard or Red-Black, depending on `mg%smoother_type`) for the variable-coefficient Poisson equation $\nabla \cdot (D \nabla \phi) = f$ on the data within a single box `id`. It updates the solution variable (`mg_iphi`) using values from neighboring cells, the source term (`mg_irhs`), and the spatially varying coefficient $D$ (from `mg_iveps`). The coefficient $D$ at cell faces is computed using harmonic averaging: $D_{face} = 2 D_0 D_{neighbor} / (D_0 + D_{neighbor})$.
  - **Arguments:** (Standard signature for `mg_box_gsrb`)

- **`box_vlpl(mg, id, nc, i_out)`:**
  - **Description:** This subroutine applies the discretized variable-coefficient Laplacian operator $\nabla \cdot (D(\mathbf{x}) \nabla \phi)$ to the solution variable `phi` (stored at index `mg_iphi`) within box `id`. The result of the operation is stored in the `i_out` component of the box's `cc` array. It also uses harmonic averaging for the coefficient $D$ at cell faces.
  - **Arguments:** (Standard signature for `mg_box_op`)

### Commented-Out Code
- **`vlpl_prolong`:** The file contains a commented-out subroutine `vlpl_prolong`. This indicates that a specialized prolongation operator, potentially designed to be more accurate for variable coefficients by taking $D(\mathbf{x})$ into account during interpolation, was developed or considered. However, it is not currently part of the active code.

## Important Variables/Constants

- **Coefficient $D(\mathbf{x})$ Storage:** The spatially varying coefficient $D$ is not a module-level variable but is expected to be provided by the user within the `mg%boxes(id)%cc(..., mg_iveps)` array for each cell.
- **Harmonic Averaging:** The discretization of the $\nabla \cdot (D \nabla \phi)$ term critically relies on harmonic averaging of $D$ values from adjacent cells to determine the effective coefficient at cell interfaces. This is a standard practice to maintain physical realism (e.g., flux continuity) when coefficients change abruptly.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, `dp`, `NDIM`, variable indices (`mg_iphi`, `mg_irhs`, `mg_iveps`), geometry types (`mg_cartesian`), smoother types (`mg_smoother_gs`, `mg_smoother_gsrb`), and boundary condition type constants (`mg_bc_neumann`).
- **External Libraries:**
  - **`cpp_macros.h`:** Used for conditional compilation based on `NDIM` to provide dimension-specific implementations of the stencil operations.
- **Interactions with Other Components:**
  - **`m_multigrid`:** This is the primary consumer. `vlaplacian_set_methods` is called (typically via `mg_set_methods` in `m_multigrid` when `mg%operator_type` is `mg_vlaplacian`) to register `box_vlpl` as `mg%box_op` and `box_gs_vlpl` as `mg%box_smoother`. The generic multigrid cycling routines then invoke these procedures.
  - **`m_ghost_cells`:** Before the operator or smoother can be applied, ghost cells for the solution (`mg_iphi`), the right-hand side (`mg_irhs` for the smoother), and importantly, the variable coefficient (`mg_iveps`) must be correctly filled. The `vlaplacian_set_methods` routine configures Neumann boundary conditions for `mg_iveps` to assist in this.
  - **Application Code:** Users wishing to solve equations like $\nabla \cdot (D(\mathbf{x}) \nabla \phi) = f$ would set `mg%operator_type = mg_vlaplacian`, populate `mg_iveps` with the coefficient $D(\mathbf{x})$ at all relevant grid levels, and then use the standard multigrid solvers.

## Usage Examples

```fortran
! Conceptual Example: Setting up the multigrid solver for a variable-coefficient
! Poisson-like problem.

use m_data_structures
use m_vlaplacian
use m_multigrid      ! For the generic mg_set_methods

implicit none

type(mg_t) :: my_problem_config
integer    :: lvl, i, box_idx, num_cells
integer    :: D_coeff_idx ! Should be mg_iveps

! 1. Initialize 'my_problem_config' (grid, MPI setup, etc.)
!    call comprehensive_mg_initialization(my_problem_config) ! Placeholder
my_problem_config%geometry_type = mg_cartesian
my_problem_config%smoother_type = mg_smoother_gsrb

! 2. Specify the need for an extra variable to store D(x).
!    The index mg_iveps (usually 5) is used for this.
D_coeff_idx = mg_iveps
my_problem_config%n_extra_vars = max(my_problem_config%n_extra_vars, D_coeff_idx - mg_num_vars + 1)

! 3. After tree building and load balancing, allocate storage.
!    call mg_allocate_storage(my_problem_config) ! This will now include space for D_coeff_idx.

! 4. Populate the variable coefficient D(x) and the source term f.
!    This must be done for all locally owned boxes on all relevant levels.
do lvl = my_problem_config%lowest_lvl, my_problem_config%highest_lvl
  num_cells = my_problem_config%box_size_lvl(lvl)
  do i = 1, size(my_problem_config%lvls(lvl)%my_ids)
    box_idx = my_problem_config%lvls(lvl)%my_ids(i)
    ! Example: Set D(x) based on spatial location, and f to some function
    ! my_problem_config%boxes(box_idx)%cc(..., D_coeff_idx) = calculate_D_at_cell(...)
    ! my_problem_config%boxes(box_idx)%cc(..., mg_irhs)    = calculate_f_at_cell(...)
  end do
end do

! 5. Set the operator type to indicate a variable-coefficient Laplacian problem.
my_problem_config%operator_type = mg_vlaplacian

! 6. Associate the vlaplacian operator and smoother with the mg_t object.
!    This call will also configure Neumann BCs for D_coeff_idx (mg_iveps).
call mg_set_methods(my_problem_config)

! Now, 'my_problem_config' is ready. When multigrid routines like mg_fas_fmg
! are invoked, they will use 'box_vlpl' for operator applications and
! 'box_gs_vlpl' for smoothing.
! It's crucial that ghost cells for mg_iphi, mg_irhs, AND D_coeff_idx (mg_iveps)
! are correctly filled during the multigrid cycles. The Neumann BC setup for
! mg_iveps by vlaplacian_set_methods aids in this process for the coefficient.
```
