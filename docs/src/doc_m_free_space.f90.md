# `m_free_space.f90`

## Overview

The `m_free_space` module provides a specialized solver for Poisson's equation ($\nabla^2 \phi = \rho$) in a three-dimensional Cartesian domain under free-space (open or unbounded) boundary conditions. This is particularly useful when the influence of charges or sources extends to infinity, and a fixed potential or zero-flux condition at a finite boundary is inappropriate.

The core idea is to use a Fast Fourier Transform (FFT) based Poisson solver (from an external `poisson_solver` module) on a chosen, potentially coarser, grid level (`fft_lvl`) of the multigrid hierarchy. This FFT solve computes the potential due to the overall source distribution as if it were in infinite space. The potential values from the faces of this FFT solution domain are then extracted and used as Dirichlet boundary conditions for the main multigrid solver, which handles the solution on finer grid levels or refines the solution on the `fft_lvl` itself.

**Note:** This module's functionality is only available when the library is compiled for `NDIM == 3`.

## Key Components

### Modules

- **`m_free_space`:** Encapsulates the logic for the free-space Poisson solver.

### Conditional Compilation
The entire module is conditionally compiled using `#if NDIM == 3 ... #endif`. If `NDIM` is not 3, this module provides no public interface or functionality.

### Derived Type Definition

- **`mg_free_bc_t` (private):** This type stores the state and data required for the free-space boundary condition calculation.
  - `initialized (logical)`: A flag, `.true.` if the FFT-derived boundary data has been computed and is ready.
  - `fft_lvl (integer)`: The multigrid level on which the FFT-based global solve was performed.
  - `rhs (real(dp), allocatable, :,:,:)`: An array holding the right-hand side (source term $\rho$) for the FFT solver, potentially padded for FFT requirements and aggregated from all MPI processes.
  - `karray (real(dp), pointer, :)`: A pointer to the Green's function kernel in Fourier space, obtained from the `poisson_solver` module.
  - `inv_dr(2, mg_num_neighbors) (real(dp))`: Inverse grid spacings for the two perpendicular directions on each of the six boundary faces of the `fft_lvl` domain. Used for interpolation.
  - `r_min(2, mg_num_neighbors) (real(dp))`: Minimum coordinates for the two perpendicular directions on each boundary face. Used for interpolation.
  - `bc_x0, bc_x1, bc_y0, bc_y1, bc_z0, bc_z1 (real(dp), allocatable, :,:)`: Arrays storing the computed potential values on the six faces (low/high x, y, z) of the `fft_lvl` domain. These serve as the Dirichlet data for the multigrid solver.

### Module Variable

- **`free_bc (type(mg_free_bc_t))` (private):** A module-level variable of type `mg_free_bc_t` that holds the persistent state of the free-space boundary calculation across calls.

### Public Subroutines

- **`mg_poisson_free_3d(mg, new_rhs, max_fft_frac, fmgcycle, max_res)`:**
  - **Description:** This is the primary public routine for solving the 3D free-space Poisson problem.
    1.  It first checks that the geometry is Cartesian and the operator is Laplacian.
    2.  It determines an optimal `fft_lvl` for the FFT solve. This level is chosen such that the number of unknowns on it is less than or equal to `max_fft_frac` times the total number of unknowns on the leaf nodes of the entire grid.
    3.  If `new_rhs` is `.true.` or if the determined `fft_lvl` has changed since the last call (or if it's the first call):
        *   The source term `mg_irhs` is restricted from finer levels down to `fft_lvl`.
        *   The Green's function kernel `free_bc%karray` is obtained by calling `createKernel` from the `poisson_solver` module.
        *   The RHS for the FFT solver is prepared in `free_bc%rhs` by gathering (MPI_Allreduce) and scaling (`rhs_fac = -1/(4\pi)`) the source term from `mg%boxes(:)%cc(..., mg_irhs)` on `fft_lvl`.
        *   The external `PSolver` routine (from `poisson_solver`) is called to compute the potential on the `fft_lvl` grid using FFT methods.
        *   The potential values on the six faces of this `fft_lvl` domain are extracted from the `PSolver` output and stored in `free_bc%bc_x0`, etc.
        *   These boundary values are registered with the multigrid framework via `mg_phi_bc_store`.
        *   The solution from `PSolver` is used as an initial guess for `mg_iphi` on `fft_lvl`, and then restricted/prolonged to populate other levels as an initial guess.
    4.  The `mg%bc(n, mg_iphi)%boundary_cond` procedure pointer is set to `ghost_cells_free_bc` for all faces, ensuring that subsequent ghost cell filling operations use the FFT-derived boundary data.
    5.  If `fft_lvl` is coarser than the finest grid level (`mg%highest_lvl`), a standard multigrid cycle (`mg_fas_fmg` if `fmgcycle` is `.true.`, else `mg_fas_vcycle`) is performed to refine the solution.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `new_rhs (logical, intent(in))`: Set to `.true.` if the right-hand side (source term in `mg_irhs`) has changed or if this is the first call. If `.false.`, previously computed boundary data might be reused if `fft_lvl` is unchanged, and only multigrid refinement cycles are performed.
    - `max_fft_frac (real(dp), intent(in))`: A fraction (0.0-1.0) controlling the maximum relative size of the FFT grid compared to the total grid.
    - `fmgcycle (logical, intent(in))`: If `.true.`, an FMG cycle is performed by the multigrid solver. If `.false.`, a V-cycle is performed.
    - `max_res (real(dp), intent(out), optional)`: If present, returns the maximum residual after the multigrid solve.

### Private Helper Subroutines

- **`ghost_cells_free_bc(box, nc, iv, nb, bc_type, bc)`:**
  - **Description:** This subroutine is assigned to the `boundary_cond` procedure pointer in `mg%bc`. When ghost cells are filled for a `box` at a physical boundary, this routine is called. It sets the `bc_type` to `mg_bc_dirichlet` and calculates the potential values for the ghost cells (`bc`) by calling `interp_bc` to interpolate from the appropriate pre-computed boundary plane (`free_bc%bc_x0`, etc.).
  - **Arguments:** Standard signature for `mg_subr_bc`.

- **`interp_bc(current_free_bc_data, nb, x1, x2)` elemental function, result(val):**
  - **Description:** Performs bilinear interpolation to find the potential `val` at a given point (`x1`, `x2` representing coordinates on a plane) on a specific boundary face `nb`. It uses the pre-calculated potential values stored in the `bc_x0`...`bc_z1` arrays within the `current_free_bc_data` (an instance of `mg_free_bc_t`).
  - **Arguments:**
    - `current_free_bc_data (type(mg_free_bc_t), intent(in))`: The stored free boundary condition data.
    - `nb (integer, intent(in))`: The neighbor index indicating which boundary face.
    - `x1, x2 (real(dp), intent(in))`: Coordinates on the boundary plane.
    - `val (real(dp))`: The interpolated potential value.

## Important Variables/Constants

- **`rhs_fac (real(dp), parameter)`:** A scaling factor of `-1.0 / (4 * acos(-1.0_dp))` (i.e., $-1/(4\pi)$) applied to the source term `mg_irhs` before passing it to the `PSolver`. This is because the Green's function for the Laplacian $\nabla^2 G = \delta(\mathbf{r})$ is typically $-1/(4\pi |\mathbf{r}|)$.
- **`free_bc%fft_lvl`**: Determines the grid level for the FFT-based global solve.
- **`mg_iphi`**: Stores the solution (potential $\phi$).
- **`mg_irhs`**: Stores the source term (charge density $\rho$).

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, parameter constants like `mg_iphi`, `mg_irhs`, `mg_num_neighbors`, `mg_neighb_dim`, `mg_bc_dirichlet`, etc.
- **External Libraries/Modules:**
  - **`mpi`:** Used in `mg_poisson_free_3d` for `MPI_ALLREDUCE` to sum the RHS across processes for the FFT solver.
  - **`poisson_solver` (external module):** This is a critical dependency. `m_free_space` calls `createKernel` to get the Green's function in Fourier space and `PSolver` to perform the actual FFT-based solution of the Poisson equation on the `fft_lvl` grid.
  - **`m_multigrid`:** Uses `mg_fas_fmg` and `mg_fas_vcycle` for the main multigrid solution cycles, and `mg_phi_bc_store` to register the boundary conditions.
  - **`m_restrict`:** Uses `mg_restrict_lvl` to transfer the source term and initial guess solution down to coarser levels.
  - **`m_prolong`:** Uses `mg_prolong` to transfer the initial guess solution from `fft_lvl` to finer levels.
  - **`m_ghost_cells`:** While the custom `ghost_cells_free_bc` provides the values, the framework calling `mg_fill_ghost_cells_lvl` (from `m_ghost_cells`) would trigger its execution.
- **Interactions with Other Components:**
  - This module provides a specialized solution strategy for 3D Poisson equations requiring open (free-space) boundary conditions.
  - It works in conjunction with the main multigrid solvers: the FFT provides a global solution component and sets boundary conditions, while the multigrid method refines this solution locally.
  - It effectively replaces standard fixed-value or zero-flux boundary conditions at the domain edge with dynamically computed Dirichlet conditions derived from the global free-space solve.
  - The accuracy of the free-space approximation depends on the chosen `fft_lvl` and the extent to which the sources are contained within this `fft_lvl` domain.

## Usage Examples

```fortran
! Conceptual example for using the free-space Poisson solver.
! This assumes the code is compiled for NDIM = 3.

#if NDIM == 3
use m_data_structures
use m_free_space
! Ensure other necessary modules like m_multigrid, m_default_settings etc. are used for setup.

implicit none

type(mg_t) :: mg_grid_config
logical    :: is_new_problem_rhs
real(dp)   :: fft_grid_size_fraction
logical    :: perform_fmg_cycle
real(dp)   :: final_residual

! --- Setup mg_grid_config ---
! (This part is complex and involves initializing the mg_t object,
!  e.g., using routines from m_default_settings, m_build_tree, m_allocate_storage.
!  It's crucial that mg_grid_config%operator_type is set to mg_laplacian,
!  and geometry is mg_cartesian.)
!
! Example:
! call mg_init_default(mg_grid_config) ! (Hypothetical setup)
! mg_grid_config%operator_type = mg_laplacian
! mg_grid_config%geometry_type = mg_cartesian
! call mg_build_rectangle(...)
! call mg_load_balance_serial(mg_grid_config) ! (Hypothetical setup)
! call mg_allocate_storage(mg_grid_config)
!
! Load the source term (e.g., charge density) into mg_grid_config%boxes(:)%cc(..., mg_irhs)
! --- End Setup ---


! Parameters for the free-space solver call
is_new_problem_rhs     = .true.  ! True if RHS changed or first solve
fft_grid_size_fraction = 0.5_dp  ! FFT grid has at most 50% of total unknowns
perform_fmg_cycle      = .true.  ! Use FMG for the multigrid part

call mg_poisson_free_3d(mg_grid_config, is_new_problem_rhs, &
                         fft_grid_size_fraction, perform_fmg_cycle, &
                         final_residual)

print *, "Free-space Poisson solution completed."
print *, "Final residual from multigrid part: ", final_residual
! The solution (potential phi) is now in mg_grid_config%boxes(:)%cc(..., mg_iphi)

! If you want to perform more V-cycles on an existing problem (RHS unchanged):
! is_new_problem_rhs = .false.
! perform_fmg_cycle  = .false. ! Perform a V-cycle
! call mg_poisson_free_3d(mg_grid_config, is_new_problem_rhs, &
!                         fft_grid_size_fraction, perform_fmg_cycle, &
!                         final_residual)
! print *, "Additional V-cycle completed. New residual: ", final_residual

#else
print *, "This example for m_free_space is only applicable for NDIM = 3."
#endif
```
