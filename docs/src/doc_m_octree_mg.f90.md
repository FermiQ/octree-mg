# `m_octree_mg.f90`

## Overview

The `m_octree_mg` module serves as the primary public interface and an "umbrella" module for the entire octree-mg library. Its main purpose is to simplify the integration and usage of the library in user applications. By using a single `use m_octree_mg` statement, developers gain access to the vast majority of the public entities—derived types, subroutines, functions, and named constants—from all the core functional modules of the library. This approach abstracts away the need for the user to individually `use` each specific `m_...` component module (like `m_data_structures`, `m_multigrid`, `m_build_tree`, etc.), thereby providing a more convenient and cleaner way to incorporate the octree-mg functionalities into their code.

## Key Components

### Modules

- **`m_octree_mg`:** This module itself. Its primary content is a series of `use` statements.

### Re-exported Modules and Entities

The `m_octree_mg` module uses the following specialized modules from the library. Due to the `public` attribute declared in `m_octree_mg` without an `only` clause, all public entities from these used modules are effectively re-exported and become available to any program unit that `use`s `m_octree_mg`.

The key modules included are:
-   **`m_data_structures`:** Provides all fundamental derived data types (e.g., `mg_t`, `mg_box_t`, `mg_lvl_t`), named constants for variable indices (e.g., `mg_iphi`, `mg_irhs`), operator types, boundary condition types, etc.
-   **`m_build_tree`:** Contains subroutines for constructing the octree or quadtree grid hierarchy (e.g., `mg_build_rectangle`).
-   **`m_load_balance`:** Includes routines for distributing the computational grid boxes among MPI processes (e.g., `mg_load_balance`).
-   **`m_ghost_cells`:** Manages the filling of ghost cell data, essential for parallel computations and applying boundary conditions (e.g., `mg_fill_ghost_cells_lvl`).
-   **`m_allocate_storage`:** Provides routines for allocating and deallocating memory for the grid data arrays (e.g., `mg_allocate_storage`, `mg_deallocate_storage`).
-   **`m_restrict`:** Contains subroutines for restriction operations (transferring data from finer to coarser grids, e.g., `mg_restrict_lvl`).
-   **`m_communication`:** Handles low-level MPI communication tasks, such as initializing MPI and transferring data buffers (e.g., `mg_comm_init`, `sort_and_transfer_buffers`).
-   **`m_prolong`:** Includes subroutines for prolongation operations (interpolating data from coarser to finer grids, e.g., `mg_prolong`).
-   **`m_multigrid`:** The core module implementing the multigrid solver algorithms like Full Approximation Scheme (FAS) V-cycles and FMG cycles (e.g., `mg_fas_vcycle`, `mg_fas_fmg`, `mg_set_methods`).
-   **`m_helmholtz`:** Implements the constant-coefficient Helmholtz operator and smoother (e.g., `helmholtz_set_lambda`).
-   **`m_vhelmholtz`:** Implements the variable-coefficient Helmholtz operator and smoother.
-   **`m_ahelmholtz`:** Implements the anisotropic-coefficient Helmholtz operator and smoother.
-   **`m_free_space`:** Provides a solver for 3D Poisson problems with free-space boundary conditions (e.g., `mg_poisson_free_3d`).

The `implicit none` statement is used, and then `public` makes all entities from the used modules available.

## Important Variables/Constants

This module does not define any new variables or constants itself. All such entities are inherited from the modules it `use`s, with `m_data_structures` being the primary source of named constants and data types.

## Dependencies and Interactions

-   **Primary Role as a Dependency Aggregator:** The `m_octree_mg` module's main function is to establish dependencies on nearly all other key modules within the octree-mg library.
-   **User Interaction:** This module is intended to be the main point of entry for users of the octree-mg library. By including `use m_octree_mg` in their application code, users can conveniently access the library's comprehensive suite of tools for multigrid simulations without needing to manage a long list of individual `use` statements for each component.

## Usage Examples

The following conceptual example illustrates how an end-user's application might leverage `m_octree_mg` to set up and solve a partial differential equation using the library's multigrid framework.

```fortran
! Example: A user's application program solving a PDE using the octree-mg library.

program my_pde_solver_application
  use m_octree_mg  ! Single 'use' statement provides access to the library's features.
  implicit none

  type(mg_t) :: solver_instance  ! mg_t type from m_data_structures, via m_octree_mg
  real(dp)   :: final_residual
  logical    :: has_initial_guess
  integer    :: ierr ! For MPI calls if any are direct

  ! 1. Initialize MPI environment and populate basic MPI info in solver_instance
  !    (mg_comm_init is from m_communication, available via m_octree_mg)
  call mg_comm_init(solver_instance)

  ! 2. Configure problem parameters (operator, geometry, etc.)
  !    (Constants like mg_laplacian are from m_data_structures)
  solver_instance%operator_type = mg_laplacian
  solver_instance%geometry_type = mg_cartesian
  ! ... (Other necessary initializations for domain size, box size, etc.)

  ! 3. Build the grid hierarchy
  !    (mg_build_rectangle is from m_build_tree)
  ! call mg_build_rectangle(solver_instance, domain_def, box_def, ...) ! Placeholder for actual args

  ! 4. Distribute grid boxes among MPI processes
  !    (mg_load_balance is from m_load_balance)
  call mg_load_balance(solver_instance)

  ! 5. Allocate memory for the cell-centered data arrays
  !    (mg_allocate_storage is from m_allocate_storage)
  call mg_allocate_storage(solver_instance)

  ! 6. Define boundary conditions and set initial/source terms
  !    (Accessing solver_instance%boxes(:)%cc and using constants like mg_iphi, mg_irhs)
  !    call setup_boundary_conditions_and_rhs(solver_instance) ! Placeholder

  ! 7. Set up the multigrid operator and smoother methods
  !    (mg_set_methods is from m_multigrid)
  call mg_set_methods(solver_instance)

  ! 8. Solve the system using a Full Multigrid (FMG) cycle
  !    (mg_fas_fmg is from m_multigrid)
  has_initial_guess = .false. ! Example: start from a zero initial solution
  call mg_fas_fmg(solver_instance, has_initial_guess, final_residual)

  if (solver_instance%my_rank == 0) then
    print *, "PDE solution complete. Final residual:", final_residual
  end if

  ! 9. Deallocate storage for grid data arrays
  !    (mg_deallocate_storage is from m_allocate_storage)
  call mg_deallocate_storage(solver_instance)

  ! 10. Finalize MPI environment (if initialized by this program/library)
  ! call MPI_Finalize(ierr) ! Direct MPI call

end program my_pde_solver_application
```
This example demonstrates that by simply `use m_octree_mg`, the application program gains access to necessary derived types (`mg_t`), subroutines from various modules (`mg_comm_init`, `mg_build_rectangle`, `mg_load_balance`, `mg_allocate_storage`, `mg_set_methods`, `mg_fas_fmg`, `mg_deallocate_storage`), and constants (`mg_laplacian`, `mg_cartesian`, `mg_iphi`, `mg_irhs`, `dp`). This greatly simplifies the user's code.
