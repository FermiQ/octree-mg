# `m_ahelmholtz.f90`

## Overview

The module `m_ahelmholtz` implements multigrid procedures for solving an anisotropic Helmholtz equation of the form: `div(D grad(phi)) - lambda*phi = f`. In this equation, `D` represents a spatially varying anisotropic diffusion coefficient (with components in each spatial direction), `lambda` is a scalar constant, `phi` is the unknown field, and `f` is the source term. This module is a core component of the `m_octree_mg` octree-based multigrid library.

## Key Components

### Modules

- **`m_ahelmholtz`:** This module encapsulates all the necessary subroutines and variables for defining and solving the anisotropic Helmholtz equation using the multigrid framework. It handles the setup of specific operator and smoother routines and manages the `lambda` parameter.

### Functions/Subroutines

- **`ahelmholtz_set_methods(mg)`:**
  - **Description:** Configures the multigrid structure (`mg`) with the appropriate function pointers for the anisotropic Helmholtz operator (`mg%box_op => box_ahelmh`) and the Gauss-Seidel smoother (`mg%box_smoother => box_gs_ahelmh`). It also sets up Neumann zero boundary conditions for the components of the anisotropic coefficient `D` (stored in `mg%cc` at indices `mg_iveps1`, `mg_iveps2`, `mg_iveps3`), as these are required in ghost cells.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.
- **`ahelmholtz_set_lambda(lambda)`:**
  - **Description:** Sets the global scalar value for `lambda` used in the Helmholtz equation. It enforces that `lambda` must be non-negative.
  - **Arguments:**
    - `lambda (real(dp), intent(in))`: The value for the Helmholtz coefficient `lambda`.
- **`box_gs_ahelmh(mg, id, nc, redblack_cntr)`:** (Private)
  - **Description:** Performs Gauss-Seidel relaxation (optionally red-black ordered) on a specified computational box (`id`) within the multigrid hierarchy. This routine is used as a smoother to reduce high-frequency errors.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `id (integer, intent(in))`: The identifier of the box to operate on.
    - `nc (integer, intent(in))`: The size of the computational box (number of cells in each direction).
    - `redblack_cntr (integer, intent(in))`: A counter used to determine the cell parity in red-black Gauss-Seidel.
- **`box_ahelmh(mg, id, nc, i_out)`:** (Private)
  - **Description:** Applies the discretized anisotropic Helmholtz operator (`L(phi) = div(D grad(phi)) - lambda*phi`) to the variable `phi` (stored at `mg_iphi` index) within a specified computational box (`id`). The result is stored in the `i_out` component of the box's cell-centered data array (`mg%boxes(id)%cc`).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `id (integer, intent(in))`: The identifier of the box to operate on.
    - `nc (integer, intent(in))`: The size of the computational box.
    - `i_out (integer, intent(in))`: The index within the `cc` array where the output of the operator will be stored.

## Important Variables/Constants

- **`ahelmholtz_lambda (real(dp), public, protected)`:** This module-level variable stores the scalar `lambda` (must be ≥ 0) for the Helmholtz equation: `div(D grad(phi)) - lambda*phi = f`. It is set by the `ahelmholtz_set_lambda` subroutine.

## Usage Examples

Usage of this module is typically indirect, through the higher-level multigrid solver routines provided by the `m_octree_mg` library. However, the setup process involves:

```fortran
! Conceptual example:
! Assume 'mg_object' is a properly initialized 'mg_t' type from m_data_structures,
! and 'desired_lambda' is a real(dp) variable holding the lambda value.

! 1. Set the lambda for the Helmholtz equation:
call ahelmholtz_set_lambda(desired_lambda)

! 2. Configure the multigrid methods for the Helmholtz solver:
! This assigns box_ahelmh to mg_object%box_op and
! box_gs_ahelmh to mg_object%box_smoother.
call ahelmholtz_set_methods(mg_object)

! 3. Perform multigrid cycles (e.g., V-cycle, FMG):
! These cycles will internally use the operator (mg_object%box_op) and
! smoother (mg_object%box_smoother) configured above.
! call mg_solve(mg_object, ...) ! (Actual solver call might vary)
```

## Dependencies and Interactions

- **Internal Dependencies:**
  - `m_data_structures`: This module heavily relies on `m_data_structures` for core data types like `mg_t` (the main multigrid structure), `dp` (double precision kind), `NDIM` (number of dimensions), and various integer constants defining indices for variables within cell-centered data arrays (e.g., `mg_iphi` for the solution field, `mg_irhs` for the right-hand side, `mg_iveps1`, `mg_iveps2`, `mg_iveps3` for the components of the anisotropic coefficient `D`). It also uses defined constants for boundary condition types (e.g., `mg_bc_neumann`).
- **External Libraries:**
  - `cpp_macros.h`: This is a C-preprocessor include file, likely containing project-specific macros used for conditional compilation (e.g., `#if NDIM == ...`). It's not an external library in the typical sense but a build system dependency.
- **Interactions with Other Components:**
  - The `m_ahelmholtz` module is integral to the `m_octree_mg` library. The `ahelmholtz_set_methods` subroutine directly modifies the `mg_t` object by assigning its function pointers (`box_op`, `box_smoother`) to the local `box_ahelmh` and `box_gs_ahelmh` subroutines, respectively. These function pointers are then invoked by the generic multigrid cycling routines (e.g., V-cycle, W-cycle, FMG) within the `m_octree_mg` framework.
  - The solver operates on `box` data structures, which are fundamental components of the octree-based grid managed by `m_octree_mg`.
  - The module requires that the anisotropic coefficients `D` are stored in the cell-centered data arrays of each box (at indices `mg_iveps1`, etc.) and sets Neumann boundary conditions for these coefficients to ensure correct behavior at box boundaries during operator application and smoothing.
```
