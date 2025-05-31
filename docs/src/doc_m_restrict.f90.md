# `m_restrict.f90`

## Overview

The `m_restrict` module is a fundamental part of the `m_octree_mg` multigrid library, dedicated to the restriction operation. In multigrid algorithms, restriction is the process of transferring data from a finer grid level to the immediate next coarser grid level. This data typically consists of either the solution itself or, more commonly in Full Approximation Scheme (FAS), the residual of the solution (the difference between the right-hand side and the operator applied to the current solution approximation). This module implements a full-weighting restriction method and manages the necessary MPI communication when data needs to be transferred between processes.

## Key Components

### Modules

- **`m_restrict`:** Contains subroutines for calculating restriction communication buffer sizes, the main restriction drivers, and the core logic for data restriction.

### Functions/Subroutines

- **`mg_restrict_buffer_size(mg, n_send, n_recv, dsize)`:**
  - **Description:** Calculates the maximum required MPI communication buffer sizes for data exchange during restriction operations. It iterates through the grid levels, determining for each fine level how many data packets (each corresponding to one coarse cell's worth of restricted data) need to be sent from local fine-grid children to their remote coarse-grid parents, and vice-versa for receiving.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure. Buffer needs are stored in `mg%comm_restrict`.
    - `n_send (integer, intent(out), array)`: Maximum number of items (coarse cells) to be sent by the current process.
    - `n_recv (integer, intent(out), array)`: Maximum number of items (coarse cells) to be received by the current process.
    - `dsize (integer, intent(out))`: Size of a single data element being transferred (effectively 1, as it's per coarse cell, which is `(mg%box_size/2)**NDIM` fine cells averaged).

- **`mg_restrict(mg, iv)`:**
  - **Description:** A high-level convenience routine that performs the restriction operation for a specified variable `iv` across all applicable grid levels. It iterates from the finest level (`mg%highest_lvl`) down to the level immediately above the coarsest (`mg%lowest_lvl+1`), calling `mg_restrict_lvl` for each transition.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `iv (integer, intent(in))`: Index of the variable (e.g., `mg_iphi`, `mg_ires`) to be restricted.

- **`mg_restrict_lvl(mg, iv, lvl)`:**
  - **Description:** This is the main driver subroutine for restricting the variable `iv` from a fine grid level `lvl` to its coarser parent level `lvl-1`.
    1.  For levels requiring communication (`lvl >= mg%first_normal_lvl`), it iterates through all locally owned fine grid boxes on `lvl`. If a fine box's parent on level `lvl-1` resides on a different MPI rank, the `restrict_set_buffer` routine is called. This routine computes the restricted value (one coarse cell's worth) from the fine box and packs it into a send buffer destined for the parent's rank.
    2.  It then calls `sort_and_transfer_buffers` (from `m_communication`) to execute the actual MPI send/receive operations.
    3.  Finally, it iterates through all locally owned parent boxes on the coarse level `lvl-1`. For each such parent, `restrict_onto` is called to populate its cell data by gathering contributions from all its children (whether local or received via MPI).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `iv (integer, intent(in))`: Index of the variable on the fine grid (`lvl`) to be restricted. The result is stored in the same variable index `iv` on the coarse grid (`lvl-1`).
    - `lvl (integer, intent(in))`: The fine grid level from which data is restricted.

### Private Helper Subroutines

- **`restrict_set_buffer(mg, id, iv)`:**
  - **Description:** Called for a fine grid box `id` (variable `iv`) if its parent `p_id` is on a different MPI rank. This routine computes the single coarse cell value that box `id` contributes to its parent `p_id`. This is done using a full weighting scheme (averaging the $2^{NDIM}$ fine cells within `id` that correspond to that specific coarse cell in `p_id`). This single averaged value is then packed into the send buffer for `p_rank`. The sorting index for the buffer is based on `p_id` and which child `id` is (its index `n` from 1 to `mg_num_children`).
- **`restrict_onto(mg, id, nc, iv)`:**
  - **Description:** Called for a coarse grid parent box `id` (on level `lvl-1`) to populate its data for variable `iv`. It iterates through its $2^{NDIM}$ child positions.
    - If a child `c_id` (on fine level `lvl`) is local (owned by the current MPI process), this routine directly accesses the child's data and performs the full weighting average of the $2^{NDIM}$ fine cells within `c_id` that correspond to the current coarse cell being computed in `id`.
    - If a child `c_id` is remote, the (already averaged by the sender via `restrict_set_buffer`) value for the corresponding coarse cell is unpacked from the MPI receive buffer.
    The routine correctly places these computed or received values into the appropriate cells of the coarse parent box `id`.

## Important Variables/Constants

- **Restriction Method:** The module implements **full weighting restriction**. In this method, the value of a coarse cell is the average of the values of all the $2^{NDIM}$ fine cells that it covers. For example, in 3D, this is `0.125 * sum(8 fine cell values)`.
- **Data Size (`dsize`):** In `mg_restrict_buffer_size`, `dsize` is set to `(mg%box_size/2)**NDIM`. This represents the number of cells in a *coarse* box if restriction were done box-by-box. However, the communication itself sends individual coarse cell values when parents are remote, so the effective `dsize` for `sort_and_transfer_buffers` in `mg_restrict_lvl` is `1` (representing one coarse cell value).

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Essential for `type(mg_t)`, `type(mg_box_t)`, `type(mg_lvl_t)`, the precision kind `dp`, `NDIM`, constants like `mg_num_children`, and utility functions such as `mg_ix_to_ichild` (used for buffer sorting index) and `mg_get_child_offset` (used to map child data to parent cells).
  - **`m_communication`:** Relies on the `sort_and_transfer_buffers` subroutine from this module to perform the actual MPI exchange of restricted data packets.
- **External Libraries:**
  - **`cpp_macros.h`:** Used for `DTIMES` (dimension-agnostic array slicing), `KJI_DO`/`CLOSE_DO` (looping macros), and `NDIM`-dependent conditional compilation of the averaging stencils.
- **Interactions with Other Components:**
  - **`m_multigrid`:** This is a primary interaction point. The `mg_restrict_lvl` routine is invoked by the `update_coarse` subroutine within the FAS V-cycle and FMG algorithms. In FAS, restriction is used to transfer both the current fine-grid solution approximation (which becomes the old solution on the coarse grid) and the fine-grid residual (which contributes to the coarse-grid right-hand side).
  - **`m_allocate_storage` (indirectly):** The `mg_restrict_buffer_size` routine is called during the library's setup phase (likely by `m_allocate_storage` or a similar top-level initialization routine that pre-calculates all buffer needs) to determine the maximum buffer sizes required for restriction data exchange.

## Usage Examples

Restriction is mostly an internal operation within the multigrid cycles, orchestrated by routines in `m_multigrid`. Users typically do not call `mg_restrict_lvl` directly.

```fortran
! Conceptual Example: How restriction is invoked within m_multigrid's update_coarse.

! In m_multigrid, within a routine like update_coarse(mg, fine_lvl):
! ...
! ! Assume 'mg' is type(mg_t)
! ! 'fine_lvl' is the current fine level index from which we are restricting.
! ! The solution on 'fine_lvl' is in 'mg_iphi'.
! ! The residual on 'fine_lvl' (f_f - L_f phi_f) has just been computed and is in 'mg_ires'.

! if (fine_lvl > mg%lowest_lvl) then
!   ! Restrict the current fine-grid solution (mg_iphi on fine_lvl)
!   ! to the coarse grid (mg_iphi on fine_lvl-1). This will serve as phi_coarse_old.
!   call mg_restrict_lvl(mg, mg_iphi, fine_lvl)
!
!   ! Restrict the fine-grid residual (mg_ires on fine_lvl)
!   ! to the coarse grid (mg_ires on fine_lvl-1).
!   call mg_restrict_lvl(mg, mg_ires, fine_lvl)
!
!   ! The coarse grid (level 'fine_lvl - 1') now has:
!   ! - its 'mg_iphi' field populated with the restricted fine-level solution.
!   ! - its 'mg_ires' field populated with the restricted fine-level residual.
!   ! These are then used to form the FAS coarse grid equation:
!   ! L_c(phi_c) = L_c(restrict(phi_f)) + restrict(f_f - L_f(phi_f))
! ...
! end if

! Direct call (less common for end-users):
use m_data_structures
use m_restrict

type(mg_t) :: mg_state
integer    :: level_to_restrict_from, variable_index

! ... (mg_state is fully initialized, data exists on level_to_restrict_from for variable_index)
! level_to_restrict_from = mg_state%highest_lvl
! variable_index = mg_iphi

! if (level_to_restrict_from > mg_state%lowest_lvl) then
!    call mg_restrict_lvl(mg_state, variable_index, level_to_restrict_from)
! end if
! Now, variable_index on level (level_to_restrict_from - 1) contains restricted data.
```
