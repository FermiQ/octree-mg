# `m_multigrid.f90`

## Overview

The `m_multigrid` module is the core engine of the `m_octree_mg` library, implementing the sophisticated Full Approximation Scheme (FAS) multigrid algorithms. Multigrid methods are highly efficient iterative techniques for solving systems of equations arising from the discretization of partial differential equations. This module orchestrates the fundamental multigrid operations:

-   **Smoothing:** Applying a few iterations of a basic iterative solver (e.g., Gauss-Seidel) to damp high-frequency errors on a given grid level.
-   **Restriction:** Transferring the problem (solution and residual) from a finer grid to a coarser grid.
-   **Prolongation:** Interpolating a correction computed on a coarser grid back to a finer grid.
-   **Coarse-Grid Solve:** Solving the problem on the coarsest grid, either directly or with more iterations of the smoother.

It provides implementations for V-cycles and Full Multigrid (FMG) cycles, which combine these operations to achieve rapid convergence. The module is designed to be generic, relying on procedure pointers within the `mg_t` data structure to call specific routines for the PDE operator, smoother, and prolongation method, which are set based on the physics of the problem (e.g., Laplacian, Helmholtz).

## Key Components

### Modules

- **`m_multigrid`:** Contains the implementations of FAS V-cycles, FMG cycles, method dispatching, and helper routines for the multigrid process.

### Functions/Subroutines

- **`mg_set_methods(mg)`:**
  - **Description:** A dispatcher routine that configures the multigrid structure `mg` with the appropriate operator-specific subroutines. Based on the value of `mg%operator_type` (e.g., `mg_laplacian`, `mg_helmholtz`), it calls the corresponding `xxx_set_methods` subroutine from the relevant operator module (e.g., `m_laplacian%laplacian_set_methods`). These specialized routines then assign the actual Fortran subroutines for operator application (`mg%box_op`), smoothing (`mg%box_smoother`), and potentially prolongation (`mg%box_prolong`) to the procedure pointers in the `mg_t` object. It also sets `mg%n_smoother_substeps` (e.g., 2 for Red-Black Gauss-Seidel).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure to be configured.

- **`mg_fas_fmg(mg, have_guess, max_res)`:**
  - **Description:** Implements the Full Multigrid (FMG) cycle using the Full Approximation Scheme (FAS). FMG provides a highly efficient way to obtain a solution by starting on the coarsest grid and progressively solving on finer grids, using the solution from the previous coarser grid as an initial guess.
    1.  If `have_guess` is `.false.`, the solution `mg_iphi` is initialized to zero on all levels.
    2.  The problem is restricted down to the `mg%lowest_lvl` by repeatedly calling `update_coarse`.
    3.  If `mg%subtract_mean` is true (for fully periodic problems), the mean of the source term `mg_irhs` is subtracted.
    4.  It then iterates from `mg%lowest_lvl` up to `mg%highest_lvl`:
        a.  The solution from the next coarser level (`lvl-1`) is prolonged to correct the current solution on `lvl` (`correct_children`).
        b.  A FAS V-cycle (`mg_fas_vcycle`) is performed on `lvl`.
    5.  Optionally returns the maximum absolute residual `max_res` on the finest level.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid structure.
    - `have_guess (logical, intent(in))`: If `.true.`, assumes `mg_iphi` contains an initial guess. If `.false.`, `mg_iphi` is initialized to zero.
    - `max_res (real(dp), intent(out), optional)`: If present, stores the maximum absolute residual on the finest level after the FMG cycle.

- **`mg_fas_vcycle(mg, highest_lvl, max_res, standalone)`:**
  - **Description:** Implements a single FAS V-cycle. A V-cycle consists of:
    1.  Pre-smoothing (relaxing) on the current grid (`smooth_boxes`).
    2.  Computing the residual ($r_f = f_f - L_f \phi_f$).
    3.  Restricting the residual and the current solution $\phi_f$ to the next coarser grid (`update_coarse`). The coarse grid problem becomes $L_c \phi_c = f_c$, where $f_c = L_c (\text{restrict}(\phi_f)) + \text{restrict}(r_f)$.
    4.  Recursively calling the V-cycle for the coarse grid problem (or solving directly if on the `mg%lowest_lvl` using more smoothing iterations).
    5.  Prolonging the computed correction from the coarse grid ($\phi_c^{\text{new}} - \phi_c^{\text{old}}$) and adding it to the fine grid solution (`correct_children`).
    6.  Post-smoothing on the current grid (`smooth_boxes`).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid structure.
    - `highest_lvl (integer, intent(in), optional)`: The finest level for this V-cycle. Defaults to `mg%highest_lvl`.
    - `max_res (real(dp), intent(out), optional)`: If present, stores the maximum absolute residual on `highest_lvl` after the V-cycle.
    - `standalone (logical, intent(in), optional)`: If `.true.` (default), timers are started/stopped and initial ghost cells are filled. If `.false.`, these are skipped (used when called from `mg_fas_fmg`).

- **`mg_apply_op(mg, i_out, op)`:**
  - **Description:** Applies the currently configured multigrid operator (pointed to by `mg%box_op`, or an optionally provided `op`) to the solution variable `mg_iphi`. This is done for all locally owned boxes across all grid levels from `mg%lowest_lvl` to `mg%highest_lvl`. The result of the operator $L(\phi)$ is stored in the variable indexed by `i_out`.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid structure.
    - `i_out (integer, intent(in))`: The index in the `cc` array where the operator result will be stored.
    - `op (procedure(mg_box_op), optional)`: An optional specific operator procedure to use instead of `mg%box_op`.

### Private Helper Subroutines (Selected)

- **`check_methods(mg)`:** Ensures `mg%box_op` and `mg%box_smoother` are associated; calls `mg_set_methods` if they are not. This is called at the beginning of `mg_fas_fmg` and `mg_fas_vcycle`.
- **`mg_add_timers(mg)`:** Initializes various timers (e.g., `timer_total_vcycle`, `timer_smoother`) used for profiling different parts of the multigrid cycles.
- **`update_coarse(mg, lvl)`:** Key FAS routine. Computes fine-grid residual $r_f = f_f - L_f \phi_f$ (using `residual_box`), restricts $\phi_f$ (to `mg_iphi` on coarse) and $r_f$ (to `mg_ires` on coarse). Then forms the coarse-grid RHS: $f_c = L_c \phi_c + r_c$. Stores the restricted $\phi_f$ into `mg_iold` on the coarse grid for later use in `correct_children`.
- **`correct_children(mg, lvl)`:** Key FAS routine. Computes the correction term on the coarse grid (`lvl`): $\delta_c = \phi_c^{\text{new}} - \phi_c^{\text{old}}$ (where $\phi_c^{\text{old}}$ is the restricted fine solution from `update_coarse`). This correction $\delta_c$ (stored in `mg_ires` on coarse grid) is then prolonged to the fine grid (`lvl+1`) and added to its solution: $\phi_f^{\text{new}} = \phi_f^{\text{old}} + P(\delta_c)$.
- **`smooth_boxes(mg, lvl, n_cycle)`:** Performs `n_cycle` smoothing iterations on all locally owned boxes at level `lvl`. Each iteration consists of `mg%n_smoother_substeps` calls to the `mg%box_smoother` procedure, with ghost cell updates (`mg_fill_ghost_cells_lvl`) after each substep.
- **`residual_box(mg, id, nc)`:** Calculates the residual $r = f - L\phi$ for a single box `id`, storing the result in `mg%boxes(id)%cc(..., mg_ires)`. $L\phi$ is computed using `mg%box_op`.
- **`subtract_mean(mg, iv, include_ghostcells)`:** For fully periodic boundary conditions, this routine subtracts the global mean of the variable `iv` from itself to ensure a unique solution. Uses `MPI_Allreduce` to compute the global sum.
- **`max_residual_lvl(mg, lvl)`:** Computes the maximum absolute value of the residual (stored in `mg_ires`) for all locally owned boxes on a given `lvl`.

## Important Variables/Constants

- **Timer Variables (Module-level):** `timer_total_vcycle`, `timer_total_fmg`, `timer_smoother`, `timer_smoother_gc`, `timer_coarse`, `timer_correct`, `timer_update_coarse`. These are handles for timers defined in `m_data_structures`.
- **Cycle Control Parameters (from `mg_t` type in `m_data_structures`):**
  - `mg%n_cycle_down`: Number of pre-smoothing steps.
  - `mg%n_cycle_up`: Number of post-smoothing steps.
  - `mg%max_coarse_cycles`: Maximum smoother iterations on the coarsest grid.
  - `mg%residual_coarse_rel`, `mg%residual_coarse_abs`: Relative and absolute tolerances for stopping the coarse grid solve.
  - `mg%subtract_mean`: Logical flag to enable mean subtraction for fully periodic problems.
  - `mg%n_smoother_substeps`: Number of sub-steps per smoother application (e.g., 2 for Red-Black Gauss-Seidel).

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Absolutely fundamental. Provides `type(mg_t)` and all its members, including procedure pointers (`box_op`, `box_smoother`, `box_prolong`), variable indices (`mg_iphi`, `mg_irhs`, `mg_iold`, `mg_ires`), and various control flags and parameters.
  - **`m_prolong`:** Provides the `mg_prolong` routine (used by `correct_children`) and the default prolongation method `mg_prolong_sparse`.
  - **`m_restrict`:** Provides the `mg_restrict_lvl` routine (used by `update_coarse`).
  - **`m_ghost_cells`:** Provides `mg_fill_ghost_cells_lvl`, which is called repeatedly within smoothing loops and after correction steps.
  - **Operator-Specific Modules (e.g., `m_laplacian`, `m_helmholtz`, `m_ahelmholtz`, etc.):** The `mg_set_methods` routine calls the `xxx_set_methods` from these modules to link the actual PDE operator and smoother implementations to the procedure pointers in `mg_t`.
- **External Libraries:**
  - **`mpi`:** Used for `MPI_ALLREDUCE` in `mg_fas_vcycle` (for global residual calculation) and in `subtract_mean` (for global sum).
  - **`cpp_macros.h`:** Used for the `DTIMES` macro for dimension-agnostic array indexing.
- **Interactions with Other Components:**
  - This module is the orchestrator of the multigrid solution process. It is typically invoked by higher-level modules or main application drivers that have already set up the specific PDE problem (e.g., `m_diffusion`, `m_poisson_fft`).
  - It requires a fully initialized `mg_t` structure, including:
    - A defined grid hierarchy (`m_build_tree`).
    - Load distribution (`m_load_balance`).
    - Allocated storage (`m_allocate_storage`).
    - The specific `operator_type` set.
    - Boundary conditions defined in `mg%bc`.
    - The source term $f$ loaded into `mg%boxes(:)%cc(..., mg_irhs)`.
  - The operator-specific modules provide the "physics" by supplying the discrete operator $L$ and smoother $S$ via procedure pointers. This module then applies them in the sequence defined by the V-cycle or FMG algorithm.

## Usage Examples

```fortran
! Conceptual Example: Solving a PDE using the FAS FMG cycle.

use m_data_structures
use m_multigrid
! Other modules would be needed for full setup:
! use m_default_settings, use m_build_tree, use m_load_balance,
! use m_allocate_storage, use m_laplacian (or other operator module)

implicit none

type(mg_t) :: my_pde_problem
logical    :: has_initial_guess
real(dp)   :: final_max_residual

! 1. Perform comprehensive setup of 'my_pde_problem':
!    - Initialize mg_t with default values.
!    - Set problem-specific parameters (e.g., my_pde_problem%operator_type = mg_laplacian).
!    - Build the grid tree (e.g., call mg_build_rectangle).
!    - Distribute load among MPI processes (e.g., call mg_load_balance).
!    - Allocate memory for grid data (e.g., call mg_allocate_storage).
!    - Define boundary conditions in my_pde_problem%bc.
!    - Populate the source term f in my_pde_problem%boxes(:)%cc(..., mg_irhs).
!    - If an initial guess for the solution is available, populate it in
!      my_pde_problem%boxes(:)%cc(..., mg_iphi) and set has_initial_guess = .true.
!    call setup_my_specific_pde_problem(my_pde_problem, has_initial_guess) ! Placeholder

! 2. Ensure operator-specific methods are set.
!    mg_fas_fmg and mg_fas_vcycle internally call check_methods, which calls
!    mg_set_methods if needed. So, an explicit call here is often for clarity
!    or if settings need to be confirmed before the main solve.
call mg_set_methods(my_pde_problem)

! 3. Solve the PDE using the Full Multigrid (FMG) algorithm.
call mg_fas_fmg(my_pde_problem, has_initial_guess, final_max_residual)

! The solution is now in my_pde_problem%boxes(:)%cc(..., mg_iphi).
if (my_pde_problem%my_rank == 0) then
    print *, "FMG Solution complete. Final maximum residual:", final_max_residual
end if

! Example of performing an additional V-cycle if needed:
! call mg_fas_vcycle(my_pde_problem, max_res=final_max_residual)
! if (my_pde_problem%my_rank == 0) then
!    print *, "Additional V-cycle. New maximum residual:", final_max_residual
! end if

! Example of applying the configured operator L to current phi, storing result in mg_ires:
! call mg_apply_op(my_pde_problem, mg_ires)
! This would compute L(phi) and place it into the mg_ires field of each box's cc array.
```
