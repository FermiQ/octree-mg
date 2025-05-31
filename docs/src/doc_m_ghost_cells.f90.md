# `m_ghost_cells.f90`

## Overview

The `m_ghost_cells` module plays a crucial role in the `m_octree_mg` library by managing the exchange and application of data in ghost cells (also known as halo or guard cells). Ghost cells are layers of virtual cells surrounding each computational box. They are filled with data from adjacent boxes—which may reside on the same MPI process or different MPI processes—or by applying physical or other boundary conditions (like those at refinement edges).

Populating these ghost cells correctly is essential before applying any finite difference stencils (e.g., for Laplacian, Helmholtz, or other operators) near the boundaries of a box, ensuring that the stencil has access to valid data from logical neighbors. This module orchestrates the data transfer, applies various types of boundary conditions, and handles data interpolation at refinement boundaries (where grid levels of different resolutions meet).

## Key Components

### Modules

- **`m_ghost_cells`:** Contains subroutines for calculating buffer sizes, packing data for communication, triggering MPI communication (via `m_communication`), and unpacking received data or applying boundary conditions to fill ghost cell regions.

### Functions/Subroutines

- **`mg_ghost_cell_buffer_size(mg, n_send, n_recv, dsize)`:**
  - **Description:** Calculates the maximum required MPI communication buffer sizes for ghost cell data exchange across all relevant grid levels. It performs a "dry run" of the ghost cell buffering logic for each level to count the number of data elements to be sent and received by each process. The output `n_send` and `n_recv` arrays are typically used by `m_allocate_storage` to allocate sufficiently large buffers in `mg%buf`.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure. Communication needs are stored in `mg%comm_ghostcell`.
    - `n_send (integer, intent(out), array)`: Maximum number of items to be sent by the current process to any other process.
    - `n_recv (integer, intent(out), array)`: Maximum number of items to be received by the current process from any other process.
    - `dsize (integer, intent(out))`: Size of a single data element being transferred (typically `nc**(NDIM-1)` for a face).

- **`mg_fill_ghost_cells(mg, iv)`:**
  - **Description:** A convenience routine that fills ghost cells for a specified variable `iv` across *all* active grid levels in the multigrid hierarchy, from `mg%lowest_lvl` to `mg%highest_lvl`. It does this by iteratively calling `mg_fill_ghost_cells_lvl`.
  - **Arguments:**
    - `mg (type(mg_t))`: The multigrid data structure.
    - `iv (integer, intent(in))`: Index of the variable (e.g., `mg_iphi`, `mg_irhs`) for which ghost cells are to be filled.

- **`mg_fill_ghost_cells_lvl(mg, lvl, iv)`:**
  - **Description:** This is the core routine for populating ghost cells for a given variable `iv` at a specific grid `lvl`. The process involves:
    1.  Iterating through all locally owned boxes (`my_ids`) on the level to prepare data that needs to be sent to neighboring boxes on other MPI processes (`buffer_ghost_cells`).
    2.  Iterating through local boxes at refinement boundaries (`my_ref_bnds` on `lvl-1`) to prepare data for finer neighbors on `lvl` that are on other MPI processes (`buffer_refinement_boundaries`).
    3.  Invoking `sort_and_transfer_buffers` (from `m_communication`) to perform the actual MPI send/receive operations.
    4.  Iterating again through local boxes on the level to fill their ghost cells (`set_ghost_cells`). This step uses the data received from MPI, or copies data directly from local neighbors, or applies physical/refinement boundary conditions.
  - **Arguments:**
    - `mg (type(mg_t))`: The multigrid data structure.
    - `lvl (integer, intent(in))`: The grid level on which to fill ghost cells.
    - `iv (integer, intent(in))`: Index of the variable to be processed.

- **`mg_phi_bc_store(mg)`:**
  - **Description:** An optimization routine that pre-calculates and stores the ghost cell values resulting from physical boundary conditions for the primary solution variable (`mg_iphi`). These values are stored into the ghost cell region normally associated with the `mg_irhs` variable within each box's `cc` array. If `mg%phi_bc_data_stored` is `.true.`, `set_ghost_cells` will use this pre-computed data when `iv == mg_iphi`, avoiding redundant boundary condition calculations.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.

### Private Helper Subroutines

The module contains numerous private helper routines. Some of the key ones include:
- `mg_phi_bc_store_lvl`: Worker routine for `mg_phi_bc_store` for a single level.
- `buffer_ghost_cells`: Identifies and packs data from a local box that needs to be sent to its off-process neighbors.
- `buffer_refinement_boundaries`: Packs data from a coarse box for its off-process fine neighbors at a refinement interface.
- `set_ghost_cells`: Orchestrates the filling of ghost cells for a single box by determining the source of data for each neighbor face (MPI buffer, local copy, physical BC, refinement BC).
- `fill_refinement_bnd`: Specifically handles filling ghost cells for faces at a refinement boundary, using either a user-defined procedure (`mg%bc(nb,iv)%refinement_bnd`) or the default `sides_rb` interpolation.
- `copy_from_nb`: Directly copies data from a neighboring box on the same MPI process.
- `buffer_for_nb` / `buffer_for_fine_nb`: Lower-level routines to place data into send buffers.
- `fill_buffered_nb`: Retrieves data from a receive buffer and places it into a box's ghost cells.
- `box_gc_for_neighbor` / `box_gc_for_fine_neighbor`: Extract data from the interior of a box that corresponds to a neighbor's ghost cell layer.
- `box_get_gc` / `box_set_gc`: Low-level accessors to get or set data in the 1-cell wide ghost layer of a `mg_box_t%cc` array (e.g., at index 0 or `nc+1`).
- `bc_to_gc`: Applies standard boundary conditions (Dirichlet, Neumann, continuous) by calculating ghost cell values based on interior points and the specified BC type and value.
- `sides_rb`: The default routine for interpolating data at refinement boundaries, designed to preserve diffusive fluxes.

## Important Variables/Constants

- **Ghost Cell Width:** The implementation implicitly assumes a ghost cell layer that is one cell wide. This is evident in routines like `box_set_gc` and `bc_to_gc` which access indices `0` and `nc+1` (where `1..nc` are interior cells).
- **Boundary Condition Types:** Constants like `mg_bc_dirichlet`, `mg_bc_neumann`, `mg_bc_continuous` (defined in `m_data_structures`) are interpreted by `bc_to_gc` to apply the correct formula.
- **`mg%phi_bc_data_stored (logical)`:** A flag in `mg_t` indicating if `mg_phi_bc_store` has been called.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** This is a fundamental dependency. The module extensively uses `type(mg_t)`, `type(mg_box_t)`, `type(mg_lvl_t)`, various named constants for neighbor indexing (e.g., `mg_num_neighbors`, `mg_neighb_rev`), variable indices (`mg_iphi`, `mg_irhs`), boundary condition types, and utility functions like `mg_get_child_offset`.
  - **`m_communication`:** Relies on `sort_and_transfer_buffers` from this module to perform the actual MPI communication of packed ghost cell data.
- **External Libraries:**
  - **`cpp_macros.h`:** Used for the `DTIMES` macro, enabling dimension-agnostic array operations in routines like `set_rhs` (though `set_rhs` is not in this module, the same macro is used here).
- **Interactions with Other Components:**
  - **Operator and Solver Modules (e.g., `m_ahelmholtz`, `m_laplacian`, `m_multigrid`):** Any module that applies numerical stencils to the grid data typically requires `mg_fill_ghost_cells_lvl` to be called beforehand to ensure that data in the ghost cell regions is valid and up-to-date. This is essential for computing correct values at cells near box boundaries.
  - **Multigrid Operations (`m_multigrid`):** Ghost cell updates are integral to the multigrid cycle:
    - Before smoother applications on each level.
    - Before calculating residuals on each level.
    - After prolongation to ensure the newly interpolated values on the fine grid can be used to correctly fill ghost cells at boundaries, especially refinement boundaries.
  - **Boundary Condition Customization:** Users can provide custom routines for physical boundary conditions and refinement boundary interpolation via procedure pointers in `mg%bc(nb, iv)%boundary_cond` and `mg%bc(nb, iv)%refinement_bnd` respectively. These are called by `set_ghost_cells` and `fill_refinement_bnd`. Modules like `m_free_space` utilize this to implement specialized boundary behaviors.
  - **Initialization and Allocation (`m_allocate_storage`):** The `mg_ghost_cell_buffer_size` routine is typically called during the setup phase (e.g., from `m_allocate_storage`) to determine the necessary sizes for communication buffers.

## Usage Examples

```fortran
! Conceptual example of filling ghost cells for a specific variable at a given level.

use m_data_structures
use m_ghost_cells
! Assume 'mg_instance' is a fully initialized 'type(mg_t)' variable.
! Assume 'current_level_idx' is the integer index of the level being processed.
! Assume 'variable_to_update' is the index (e.g., mg_iphi) of the variable
! in mg_instance%boxes(:)%cc whose ghost cells need updating.

implicit none
type(mg_t) :: mg_instance
integer    :: current_level_idx
integer    :: variable_to_update

! ... (mg_instance is set up, tree created, storage allocated, data initialized) ...

current_level_idx = mg_instance%highest_lvl ! Example: highest level
variable_to_update = mg_iphi             ! Example: the solution variable

! Fill ghost cells for 'variable_to_update' at 'current_level_idx'
call mg_fill_ghost_cells_lvl(mg_instance, current_level_idx, variable_to_update)

! At this point, numerical operators (smoothers, residual calculations) can be applied
! to the boxes on 'current_level_idx', as their ghost cells are now consistent
! with neighboring data and boundary conditions.

! --- Example of using pre-stored boundary conditions for mg_iphi ---
! This is an optimization if boundary conditions for mg_iphi are static.

! During initialization phase, after defining mg_instance%bc:
if (.not. mg_instance%phi_bc_data_stored) then
  call mg_phi_bc_store(mg_instance)
end if

! Later, when needing to fill ghost cells for mg_iphi:
call mg_fill_ghost_cells_lvl(mg_instance, current_level_idx, mg_iphi)
! The call to set_ghost_cells inside mg_fill_ghost_cells_lvl will detect that
! mg_phi_bc_data_stored is true and iv is mg_iphi, and will use the pre-calculated
! values from the mg_irhs ghost region for physical boundaries.
```
