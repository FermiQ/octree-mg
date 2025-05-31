# `m_allocate_storage.f90`

## Overview

The `m_allocate_storage` module serves as the primary memory management unit within the `m_octree_mg` (Octree Multigrid) library. It is responsible for the dynamic allocation and deallocation of memory for essential data structures, particularly for the cell-centered data within computational grid boxes (`mg%boxes(id)%cc`), arrays describing the multigrid level properties (`mg%lvls`), and the communication buffers (`mg%buf`) required for parallel processing tasks like ghost cell exchange, restriction, and prolongation.

## Key Components

### Modules

- **`m_allocate_storage`:** This module provides the public interface for allocating and freeing memory associated with a multigrid (`mg_t`) object.

### Functions/Subroutines

- **`mg_allocate_storage(mg)`:**
  - **Description:** This subroutine allocates memory for the primary data arrays once the multigrid tree structure (`mg%tree_created`) has been established. Specifically, it allocates:
    1.  The cell-centered data array `cc` for each locally owned computational box (`mg%boxes(id)`) on each multigrid level. The dimensions of `cc` depend on the box size (`nc`), the number of dimensions (`NDIM` via `DTIMES` macro), and the total number of variables (`mg_num_vars + mg%n_extra_vars`). Initial values are set to zero.
    2.  Communication buffers (`mg%buf`) for each potential MPI rank. The size of these buffers is determined by querying routines from `m_ghost_cells`, `m_restrict`, and `m_prolong` modules to ensure they are large enough for the respective communication patterns.
    It sets `mg%is_allocated` to true upon successful completion.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure for which storage is to be allocated.
- **`mg_deallocate_storage(mg)`:**
  - **Description:** This subroutine deallocates all dynamically allocated memory associated with the provided multigrid structure (`mg`). This includes:
    1.  The `cc` data arrays within all `mg%boxes`.
    2.  The send/receive buffers and index arrays within `mg%buf`.
    3.  Various arrays associated with each multigrid level stored in `mg%lvls(lvl)%...` (e.g., `ids`, `leaves`, `parents`, `my_ids`).
    4.  Communication partner lists for restriction, prolongation, and ghost cells.
    It also resets `mg%is_allocated` to `false` and `mg%n_boxes` to 0.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure whose storage is to be deallocated.

## Important Variables/Constants

This module primarily acts on variables defined within the `mg_t` type (from `m_data_structures`). Key data members from `mg_t` that influence allocation include:

- **`mg%box_size_lvl(lvl)`:** An array determining the number of cells (`nc`) per side for boxes at a given multigrid level `lvl`.
- **`mg_num_vars`:** A constant (from `m_data_structures`) specifying the base number of physical variables stored per cell.
- **`mg%n_extra_vars`:** A variable within `mg_t` allowing for additional user-defined variables per cell.
- **`mg%n_cpu`:** The number of MPI processes, used for sizing communication buffer arrays.
- **`mg%lowest_lvl`, `mg%highest_lvl`:** Define the range of multigrid levels to allocate.
- **`mg%lvls(lvl)%my_ids`:** An array of box IDs owned by the current process at a given level, for which `cc` data needs to be allocated.

The actual sizes of communication buffers are determined by functions imported from `m_ghost_cells`, `m_restrict`, and `m_prolong`.

## Usage Examples

The subroutines in this module are typically called after the basic multigrid hierarchy and box distribution have been determined, and before numerical operations begin. Conversely, deallocation occurs when the multigrid structure is no longer needed.

```fortran
! Conceptual example of allocation and deallocation:

! Assume 'my_mg_structure' is a variable of type(mg_t)
! and its tree (box hierarchy, levels, etc.) has been created.
! For example, after 'call mg_create_tree(my_mg_structure, ...)'

! Before performing operations, allocate storage:
if (my_mg_structure%tree_created .and. .not. my_mg_structure%is_allocated) then
    call mg_allocate_storage(my_mg_structure)
    print *, "Storage allocated."
else
    if (.not. my_mg_structure%tree_created) print *, "Tree not created yet!"
    if (my_mg_structure%is_allocated) print *, "Storage already allocated!"
end if

! ...
! Perform multigrid computations using the allocated mg%boxes(id)%cc arrays
! and communication buffers mg%buf.
! ...

! When finished, deallocate the storage:
if (my_mg_structure%is_allocated) then
    call mg_deallocate_storage(my_mg_structure)
    print *, "Storage deallocated."
else
    print *, "Storage was not allocated or already deallocated."
end if

```

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** This is a fundamental dependency. `m_allocate_storage` extensively uses and manipulates the `mg_t` derived type and its members (e.g., `boxes`, `lvls`, `buf`, `is_allocated`, `tree_created`, `box_size_lvl`, `mg_num_vars`, `n_extra_vars`, `phi_bc_data_stored`, `n_boxes`). The `DTIMES` macro, likely defined to handle array dimensioning based on `NDIM`, also originates from or is used in conjunction with this module.
  - **`m_ghost_cells` (via `mg_ghost_cell_buffer_size`):** Called by `mg_allocate_storage` to determine the required size for ghost cell communication buffers.
  - **`m_restrict` (via `mg_restrict_buffer_size`):** Called by `mg_allocate_storage` to determine the required size for restriction operation communication buffers.
  - **`m_prolong` (via `mg_prolong_buffer_size`):** Called by `mg_allocate_storage` to determine the required size for prolongation operation communication buffers.
- **External Libraries:**
  - **`cpp_macros.h`:** This preprocessor include is used, likely for the `DTIMES` macro which appears to be used for setting array bounds based on the number of spatial dimensions (`NDIM`).
- **Interactions with Other Components:**
  - **Tree Construction Modules (e.g., `m_tree_build` - assumed):** Modules responsible for creating the multigrid hierarchy (`mg%tree_created = .true.`) must be called before `mg_allocate_storage`.
  - **Parallel Communication Modules (e.g., `m_mpi_common`, `m_ghost_cells`, `m_restrict`, `m_prolong`):** These modules will utilize the communication buffers (`mg%buf`) allocated by `mg_allocate_storage`.
  - **Solver and Operator Modules (e.g., `m_ahelmholtz`, `m_poisson`):** These modules operate on the `mg%boxes(id)%cc` data arrays that are allocated by `mg_allocate_storage`.
  - **Initialization/Finalization routines:** The main program or higher-level control modules would invoke `mg_allocate_storage` during setup and `mg_deallocate_storage` during cleanup.
```
