# `m_prolong.f90`

## Overview

The `m_prolong` module is a core component of the `m_octree_mg` multigrid library, responsible for the prolongation operation, also commonly known as interpolation. In the context of multigrid methods, prolongation is the process of transferring data (typically a correction computed on a coarser grid) to a finer grid. This interpolated correction is then used to update the solution on the finer grid. This module provides the necessary routines to perform this interpolation, including handling the communication required when the coarse grid box and its corresponding fine grid children reside on different MPI processes.

## Key Components

### Modules

- **`m_prolong`:** Contains subroutines for calculating prolongation communication buffer sizes, the main prolongation driver, and specific interpolation methods.

### Functions/Subroutines

- **`mg_prolong_buffer_size(mg, n_send, n_recv, dsize)`:**
  - **Description:** Calculates the required MPI communication buffer sizes for data exchange during prolongation operations. It appears to derive these sizes from the communication patterns already established for restriction operations (`mg%comm_restrict`), as prolongation often involves communication with the same set of processes but in the reverse direction (coarse-to-fine vs. fine-to-coarse). The `dsize` is based on the number of cells in a fine grid box, as fine grid data (after interpolation on the coarse data owner's side) is sent.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure. Buffer needs are stored in `mg%comm_prolong`.
    - `n_send (integer, intent(out), array)`: Maximum number of items to be sent by the current process.
    - `n_recv (integer, intent(out), array)`: Maximum number of items to be received by the current process.
    - `dsize (integer, intent(out))`: Size of a single data element being transferred (number of cells in a fine box: `mg%box_size**NDIM`).

- **`mg_prolong(mg, lvl, iv, iv_to, method, add)`:**
  - **Description:** This is the main public routine for performing prolongation. It transfers data from variable `iv` on the coarse grid level `lvl` to variable `iv_to` on the fine grid level `lvl+1`.
    1.  For levels requiring communication (`lvl >= mg%first_normal_lvl-1`), it iterates through locally owned coarse boxes. If a child box on the fine level belongs to a different MPI rank, `prolong_set_buffer` is called to perform the interpolation (using `method`) and pack the result into a send buffer.
    2.  It then calls `sort_and_transfer_buffers` (from `m_communication`) to execute the MPI data exchange.
    3.  Finally, it iterates through locally owned fine boxes on level `lvl+1`. `prolong_onto` is called for each, which either uses locally available parent data or data from the receive buffer to perform the interpolation (via `method`) and update the `iv_to` variable.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `lvl (integer, intent(in))`: The coarse grid level to prolong from.
    - `iv (integer, intent(in))`: Index of the source variable on the coarse grid.
    - `iv_to (integer, intent(in))`: Index of the target variable on the fine grid.
    - `method (procedure(mg_box_prolong))`: The specific prolongation/interpolation procedure to be used (e.g., `mg_prolong_sparse`).
    - `add (logical, intent(in))`: If `.true.`, the prolonged value is added to the existing value in `iv_to` on the fine grid. If `.false.`, it overwrites it. (Adding is typical for FAS corrections).

- **`mg_prolong_sparse(mg, p_id, dix, nc, iv, fine)`:**
  - **Description:** A specific public prolongation method that implements multilinear interpolation (linear in 1D, bilinear in 2D, trilinear in 3D). It interpolates data from a coarse parent box `p_id` (variable `iv`) to a single fine child box, whose position relative to the parent is given by the logical offset `dix` (identifying which of the $2^{NDIM}$ children it is). The interpolated values for the child box (of size `nc`) are returned in the `fine` array. This is often the default prolongation method assigned to `mg%box_prolong`.
  - **Arguments:** (Standard signature for `mg_box_prolong` procedure pointer)
    - `mg (type(mg_t), intent(inout))`: Multigrid structure.
    - `p_id (integer, intent(in))`: ID of the parent (coarse) box.
    - `dix(NDIM) (integer, intent(in))`: Logical index offset of the child within the parent's grid structure (e.g., `[0,0]` for the first child in 2D, assuming 0-based indexing for offsets).
    - `nc (integer, intent(in))`: Size of the child (fine) box.
    - `iv (integer, intent(in))`: Index of the variable to prolong from the coarse box.
    - `fine (real(dp), intent(out), array)`: Output array to store the interpolated values for the fine child box. Its dimensions depend on `nc` and `NDIM`.

### Private Helper Subroutines

- **`prolong_set_buffer(mg, id, nc, iv, method)`:**
  - **Description:** If any children of the coarse box `id` reside on different MPI ranks, this routine calls the specified `method` (e.g., `mg_prolong_sparse`) to compute the interpolated data for those remote children. The resulting fine grid data is then packed into the appropriate send buffer in `mg%buf(child_rank)%send`.
- **`prolong_onto(mg, id, nc, iv, iv_to, add, method)`:**
  - **Description:** For a given fine grid box `id`, this routine orchestrates the actual interpolation. It determines if its parent box `p_id` is local or remote. If local, it calls `method` directly. If remote, it unpacks the previously received interpolated data from `mg%buf(parent_rank)%recv`. The result is then either added to or assigned to `mg%boxes(id)%cc(..., iv_to)`.

## Important Variables/Constants

- **Interpolation Stencils:** The `mg_prolong_sparse` routine uses hardcoded weights for multilinear interpolation. For example, in 1D, it uses weights of `0.75` and `0.25` (effectively $(c_i + c_{i-1})/2$ and $(c_i + c_{i+1})/2$ for cell centers if $c_i$ is a cell center, but adapted for cell-to-cell interpolation). Similar fixed coefficient patterns are used for 2D (bilinear) and 3D (trilinear) interpolation.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, `type(mg_lvl_t)`, the precision kind `dp`, `NDIM`, constants like `mg_num_children`, and utility functions like `mg_get_child_offset`.
  - **`m_communication`:** Relies on the `sort_and_transfer_buffers` subroutine from this module to perform the actual MPI exchange of interpolated fine-grid data when needed.
- **External Libraries:**
  - **`cpp_macros.h`:** Used for `DTIMES` (dimension-agnostic array slicing) and `NDIM`-dependent conditional compilation of interpolation stencils.
- **Interactions with Other Components:**
  - **`m_multigrid`:** This is a primary interaction. The `mg_prolong` routine (or more accurately, a specific method like `mg_prolong_sparse` that is assigned to the `mg%box_prolong` procedure pointer) is invoked by the `correct_children` subroutine within the FAS V-cycle and FMG algorithms. This step applies the coarse-grid correction to the fine grid.
  - **`m_allocate_storage` (indirectly):** The `mg_prolong_buffer_size` routine is called during the library's setup phase (likely by `m_allocate_storage` or a similar top-level initialization routine that pre-calculates all buffer needs) to determine the maximum buffer sizes required for prolongation data exchange.
  - **`m_ghost_cells`:** After data is prolonged to a fine level and corrections are applied, the ghost cells for that fine level typically need to be updated (using routines from `m_ghost_cells`) before subsequent operations like smoothing can be performed correctly.

## Usage Examples

Prolongation is mostly an internal operation within the multigrid cycles. Users typically don't call `mg_prolong` directly but rely on it being called by `mg_fas_vcycle` or `mg_fas_fmg`.

```fortran
! Conceptual Example: How prolongation is invoked within m_multigrid's correct_children.

! In m_multigrid, within a routine like correct_children(mg, coarse_lvl):
! ...
! ! Assume 'mg' is type(mg_t)
! ! 'coarse_lvl' is the current coarse level index.
! ! The correction term (phi_coarse_new - phi_coarse_old) has been computed and
! ! is stored in the 'mg_ires' field of the coarse grid boxes on 'coarse_lvl'.
!
! ! The goal is to prolong this correction from 'mg_ires' on 'coarse_lvl'
! ! and add it to 'mg_iphi' on the fine level ('coarse_lvl + 1').
!
! if (coarse_lvl < mg%highest_lvl) then
!   ! mg%box_prolong would point to a specific method like mg_prolong_sparse
!   call mg_prolong(mg, coarse_lvl, mg_ires, mg_iphi, &
!                   mg%box_prolong, add=.true.)
!
!   ! Following this, m_multigrid would typically fill ghost cells for mg_iphi
!   ! on level 'coarse_lvl + 1'.
!   ! call mg_fill_ghost_cells_lvl(mg, coarse_lvl + 1, mg_iphi)
! end if
! ...

! The default prolongation method mg_prolong_sparse can also be called directly
! if needed for specific purposes, though this is less common for end-users:
use m_data_structures
use m_prolong

type(mg_t) :: mg_instance
integer    :: parent_box_id, child_offset(1), fine_box_size, variable_index
real(dp), allocatable :: interpolated_data_for_child(:)

! ... (mg_instance, parent_box_id, child_offset, fine_box_size, variable_index are set up)
! fine_box_size = mg_instance%box_size_lvl(mg_instance%boxes(parent_box_id)%lvl + 1)
! allocate(interpolated_data_for_child(fine_box_size)) ! Assuming 1D for simplicity

! call mg_prolong_sparse(mg_instance, parent_box_id, child_offset, &
!                        fine_box_size, variable_index, interpolated_data_for_child)

! 'interpolated_data_for_child' now holds the data for one fine child.
! if (allocated(interpolated_data_for_child)) deallocate(interpolated_data_for_child)

```
