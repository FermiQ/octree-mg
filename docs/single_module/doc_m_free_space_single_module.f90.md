# `m_free_space.f90` (from `single_module` directory)

## Overview

The `m_free_space` module, specifically the version found in the `single_module` directory, provides specialized functionality to solve 3D Poisson's equation ($\nabla^2 \phi = \rho$) under free-space (open or unbounded) boundary conditions. It is designed to be compiled and used in conjunction with the `m_octree_mg_3d` module, which is the single-file, 3D version of the core `octree-mg` library.

This module employs a hybrid solution strategy. It utilizes an external FFT-based Poisson solver (from the `Poisson_Solver` package located in `poisson_3d_fft/`) to compute the global, far-field solution on a selected (potentially coarser) grid level within the `octree-mg` hierarchy. The potential values derived from this global FFT solution are then used to set Dirichlet boundary conditions for the main `octree-mg_3d` multigrid solver. This approach allows the `octree-mg` library to handle problems where the computational domain is effectively infinite, rather than being confined by simple periodic or fixed-value boundaries. This module is only functional for 3D problems.

## Key Components

### Modules

-   **`m_free_space`**: The Fortran module defined in this file.

### Internal Derived Type

-   **`mg_free_bc_t` (private)**: This derived type is used internally to store the state and data associated with the free-space boundary condition calculation. Its key components include:
    -   `initialized (logical)`: A flag indicating whether the FFT-derived boundary data has been computed.
    -   `fft_lvl (integer)`: The level within the `m_octree_mg_3d` grid hierarchy on which the global FFT solve was performed.
    -   `rhs (real(dp), allocatable, :,:,:)`: Stores the right-hand side (source term $\rho$), restricted from the main grid and padded as necessary for the FFT solver.
    -   `karray (real(dp), pointer, :)`: A pointer to the reciprocal space Poisson kernel ($G(k)$) computed by and used by the `Poisson_Solver` module.
    -   `inv_dr(2, mg_num_neighbors) (real(dp))`, `r_min(2, mg_num_neighbors) (real(dp))`: Inverse grid spacings and minimum coordinates for the boundary planes of the FFT domain, used for interpolating boundary values.
    -   `bc_x0, bc_x1, bc_y0, bc_y1, bc_z0, bc_z1 (real(dp), allocatable, :,:)`: Arrays that store the potential values on the six faces of the FFT domain, derived from the `Poisson_Solver` solution. These values are then used as Dirichlet boundary conditions for the `m_octree_mg_3d` solver.

### Module Variable (State)

-   **`free_bc (type(mg_free_bc_t))` (private)**: A module-level variable of type `mg_free_bc_t` that holds the persistent state of the free-space boundary calculation across calls to `mg_poisson_free_3d`.

### Public Subroutines

-   **`mg_poisson_free_3d(mg, new_rhs, max_fft_frac, fmgcycle, max_res)`**:
    -   **Description:** This is the main public subroutine for solving the 3D Poisson equation with free-space boundary conditions.
        1.  **Pre-checks:** Verifies that the geometry is Cartesian (`mg%geometry_type == mg_cartesian`) and the operator is Laplacian (`mg%operator_type == mg_laplacian`).
        2.  **FFT Level Selection:** Determines an appropriate `fft_lvl` for the global FFT solve. This level is chosen such that its number of unknowns is a certain fraction (`max_fft_frac`) of the total unknowns on the finest leaf level of the `mg` structure.
        3.  **Kernel and RHS Preparation (if `new_rhs` or `fft_lvl` changed):**
            -   Restricts the source term from `mg%boxes(:)%cc(..., mg_irhs)` down to the selected `fft_lvl`.
            -   Calls `createKernel` from the `poisson_solver` module to generate the free-space Green's function in reciprocal space (`free_bc%karray`).
            -   Prepares the right-hand side `free_bc%rhs` for `PSolver` by gathering data from all MPI processes (using `MPI_ALLREDUCE`) and applying a scaling factor (`rhs_fac`).
        4.  **Global Solve:** Calls `PSolver` from the `poisson_solver` module to obtain the potential on the `fft_lvl` grid.
        5.  **Boundary Extraction:** Extracts potential values from the faces of the `fft_lvl` solution domain and stores them in `free_bc%bc_x0, ...` after interpolation.
        6.  **Boundary Condition Setup:** Calls `mg_phi_bc_store` (from `m_octree_mg_3d`) to register these computed boundary values. It also sets the `mg%bc(n, mg_iphi)%boundary_cond` procedure pointer to the local `ghost_cells_free_bc` routine.
        7.  **Initial Guess:** Uses the `PSolver` solution on `fft_lvl` as an initial guess for `mg_iphi`, then restricts and prolongs this guess to other levels of the `mg` hierarchy.
        8.  **Multigrid Solve:** If the `fft_lvl` is not the highest (finest) level in the `mg` structure, it calls the main multigrid solver (`mg_fas_fmg` or `mg_fas_vcycle` from `m_octree_mg_3d`) to refine the solution using the `octree-mg` framework with the newly defined boundary conditions.
    -   **Key Arguments:**
        -   `mg (type(mg_t), intent(inout))`: The main multigrid data structure from `m_octree_mg_3d`.
        -   `new_rhs (logical, intent(in))`: Flag indicating if the source term (`mg_irhs`) has changed or if it's the first call, triggering recalculation of the FFT solution.
        -   `max_fft_frac (real(dp), intent(in))`: Fraction (0.0-1.0) determining the relative size of the coarse grid used for the FFT solve compared to the full grid.
        -   `fmgcycle (logical, intent(in))`: If `.true.`, an FMG cycle is used for the `octree-mg` solve; otherwise, a V-cycle is used.
        -   `max_res (real(dp), intent(out), optional)`: If present, returns the maximum residual from the `octree-mg` solver.

### Private Helper Subroutines

-   **`ghost_cells_free_bc(box, nc, iv, nb, bc_type, bc)`**:
    -   **Description:** This routine is assigned to the `mg%bc(n, mg_iphi)%boundary_cond` procedure pointer. When the `m_octree_mg_3d` library fills ghost cells for a box at what it considers a physical boundary, this subroutine is called. It sets `bc_type = mg_bc_dirichlet` and uses `interp_bc` to determine the ghost cell values by interpolating from the appropriate boundary plane data (e.g., `free_bc%bc_x0`) that was pre-computed by `PSolver`.
    -   **Arguments:** Standard signature for `mg_subr_bc` (from `m_data_structures`).

-   **`interp_bc(bc_data, nb, x1, x2) elemental function result(val)`**:
    -   **Description:** An elemental function that performs bilinear interpolation. Given coordinates `x1, x2` on a specific boundary face `nb`, it interpolates the potential value `val` from the corresponding pre-calculated boundary data array (e.g., `bc_data%bc_x0`) stored in the `mg_free_bc_t` type variable.
    -   **Arguments:**
        -   `bc_data (type(mg_free_bc_t), intent(in))`: The stored free boundary condition data.
        -   `nb (integer, intent(in))`: Neighbor index indicating the boundary face.
        -   `x1, x2 (real(dp), intent(in))`: Coordinates on the boundary plane.
        -   `val (real(dp))`: The interpolated potential value.

## Important Variables/Constants

-   **`free_bc (type(mg_free_bc_t))`**: Internal module variable holding the state for the free-space solver.
-   **`rhs_fac (real(dp), parameter)`**: Value of `-1.0_dp / (4 * acos(-1.0_dp))` (i.e., $-1/(4\pi)$). This scales the source term `mg_irhs` before it's passed to the `PSolver`, aligning with the Green's function definition $\nabla^2 G = -4\pi \delta(\mathbf{r})$ implicitly used by `PSolver`.
-   **Hardcoded parameters for `PSolver` call:**
    -   `geocode = 'F'` (Free boundary conditions).
    -   `datacode = 'G'` (Global data layout for the `PSolver` call's perspective on `free_bc%rhs`).
    -   `itype_scf = 8` (Order of interpolating scaling functions for `createKernel`).
    -   `ixc = 0` (No exchange-correlation calculation within this context).

## Dependencies and Interactions

-   **`use m_octree_mg_3d`**: This is the primary dependency. It provides the `mg_t` derived type, core multigrid routines (`mg_fas_fmg`, `mg_fas_vcycle`, `mg_restrict_lvl`, `mg_prolong`, `mg_fill_ghost_cells_lvl`, `mg_phi_bc_store`, `mg_number_of_unknowns`, `mg_highest_uniform_lvl`, `mg_get_face_coords`), and all associated constants (like `mg_iphi`, `mg_irhs`, `mg_laplacian`, `mg_cartesian`, `mg_bc_dirichlet`, `mg_num_neighbors`, `mg_neighb_dim`, `mg_neighb_lowx`, etc.).
-   **`use poisson_solver`**: This module (located in the `poisson_3d_fft/` directory) is crucial. `mg_poisson_free_3d` calls `createKernel` and `PSolver` from it to perform the underlying FFT-based Poisson solve on the selected coarse `fft_lvl`.
-   **`use mpi`**: Used for MPI communication, specifically `MPI_ALLREDUCE` to sum the RHS from all processes before the `PSolver` call.
-   **Interaction with `mg_t` structure (from `m_octree_mg_3d`):**
    -   Reads the grid hierarchy, domain parameters, and the source term from `mg%boxes(:)%cc(..., mg_irhs)`.
    -   Modifies `mg%bc` by assigning `ghost_cells_free_bc` as the boundary condition handler for `mg_iphi`.
    -   Invokes core multigrid solvers (`mg_fas_fmg`, `mg_fas_vcycle`) which operate on the `mg_t` structure.
    -   The final solution is stored in `mg%boxes(:)%cc(..., mg_iphi)`.

## Note on Compilation

To use this `m_free_space` module (from the `single_module` directory), the following components must be compiled and linked:
1.  The `m_octree_mg_3d` module. This is typically obtained by compiling `single_module/m_octree_mg.f90` with the preprocessor flag `-DNDIM=3`.
2.  The `single_module/m_free_space.f90` file itself (which defines the `module m_free_space`).
3.  The `poisson_solver` module. This requires compiling `poisson_3d_fft/poisson_solver.f90` along with all the source files it `include`s (e.g., `psolver_main.f90`, `build_kernel.f90`, `fft3d.f90`, etc.).
4.  An MPI library.
5.  (If XC support were enabled in `PSolver`, though it's hardcoded off here: ABINIT's XC library and its dependencies).

## Usage Examples

```fortran
! Conceptual usage of the m_free_space module with m_octree_mg_3d

program test_free_space_solver_single_module
  use m_octree_mg_3d  ! The single-file 3D octree-mg library
  use m_free_space    ! This module for free-space boundary conditions
  use mpi             ! For MPI_Init, MPI_Finalize, etc.

  implicit none

  type(mg_t) :: mg_configuration
  logical    :: is_new_rhs
  real(dp)   :: fft_coarsening_factor
  logical    :: use_fmg_solver
  real(dp)   :: final_solver_residual
  integer    :: mpi_error

  ! 1. Initialize MPI Environment
  call MPI_Init(mpi_error)
  call mg_comm_init(mg_configuration) ! Initialize MPI info within mg_t

  ! 2. Configure mg_configuration structure
  !    (Set domain size, box size, operator_type=mg_laplacian, geometry_type=mg_cartesian etc.)
  !    call setup_mg_grid_and_parameters(mg_configuration) ! Placeholder for detailed setup
  mg_configuration%operator_type = mg_laplacian
  mg_configuration%geometry_type = mg_cartesian

  ! 3. Build the octree grid structure
  !    call mg_build_rectangle(mg_configuration, ...) ! Placeholder

  ! 4. Load balance the grid
  !    call mg_load_balance(mg_configuration) ! Placeholder

  ! 5. Allocate storage for grid data
  !    call mg_allocate_storage(mg_configuration) ! Placeholder

  ! 6. Define the source term (e.g., charge density)
  !    This involves populating mg_configuration%boxes(:)%cc(..., mg_irhs)
  !    call set_source_term(mg_configuration) ! Placeholder

  ! 7. Set parameters for the free-space solver
  is_new_rhs = .true.              ! True for the first solve or if RHS changes
  fft_coarsening_factor = 0.1_dp ! FFT grid uses <= 10% of total unknowns
  use_fmg_solver = .true.          ! Use FMG for the octree-mg part

  ! 8. Call the free-space Poisson solver
  !    This routine will internally use poisson_solver for the coarse grid
  !    and then call m_octree_mg_3d routines (like mg_fas_fmg) for finer levels.
  call mg_poisson_free_3d(mg_configuration, is_new_rhs, fft_coarsening_factor, &
                           use_fmg_solver, max_res=final_solver_residual)

  ! 9. The solution (potential phi) is now in mg_configuration%boxes(:)%cc(..., mg_iphi)
  if (mg_configuration%my_rank == 0) then
    print *, "Free-space Poisson solve complete."
    print *, "Final residual from octree-mg solver:", final_solver_residual
  end if

  ! 10. Deallocate storage
  !    call mg_deallocate_storage(mg_configuration) ! Placeholder

  ! 11. Finalize MPI
  call MPI_Finalize(mpi_error)

end program test_free_space_solver_single_module
```
