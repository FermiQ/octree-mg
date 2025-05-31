# `m_load_balance.f90`

## Overview

The `m_load_balance` module in the `m_octree_mg` library is dedicated to distributing the computational workload—represented by the individual grid boxes of the octree/quadtree structure—among the available MPI processes. The goal of load balancing is to ensure that each process receives a roughly equal share of the computational effort, which is critical for achieving good parallel efficiency and scalability. This module does not transfer the actual box data (`cc` arrays) but rather assigns an owner MPI rank (`mg%boxes(id)%rank`) to each box `id`. Subsequent operations, like memory allocation (`m_allocate_storage`) and computations, then use this rank information to operate only on locally owned boxes. The overall tree structure (connectivity) is assumed to be known by all processes.

## Key Components

### Modules

- **`m_load_balance`:** Provides different strategies for assigning MPI ranks to grid boxes.

### Functions/Subroutines

- **`mg_load_balance_simple(mg)`:**
  - **Description:** Implements a straightforward load balancing scheme.
    1.  Boxes on very coarse levels (at or below `single_cpu_lvl`, which is derived from `mg%first_normal_lvl` and `mg%lowest_lvl` where box sizes might not be uniform for standard communication) are assigned to rank 0.
    2.  For finer levels, it distributes contiguous blocks of boxes on each level as evenly as possible across the MPI ranks. The distribution relies on the Morton-like ordering of box IDs within `mg%lvls(lvl)%ids` that results from the tree construction process.
    This method is noted to be most effective for grids that are largely uniform. After assigning ranks, it calls `update_lvl_info` for each level.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure whose `mg%boxes(:)%rank` fields will be updated.

- **`mg_load_balance(mg)`:**
  - **Description:** A more refined load balancing strategy that attempts to keep parent boxes on the same rank as the majority of their children.
    1.  It iterates from the finest level (`mg%highest_lvl`) down to `single_cpu_lvl + 1`.
    2.  On each level, it first assigns ranks to parent boxes: a parent is assigned to the rank that owns the most of its children (using `most_popular` to break ties by favoring ranks with less current work).
    3.  Then, it distributes the leaf boxes on that level to balance the remaining workload among MPI ranks, considering the work already assigned via parent boxes.
    4.  Boxes on levels at or below `single_cpu_lvl` are assigned to a single `coarse_rank`, which is determined by the most prevalent rank on the `single_cpu_lvl + 1` level.
    After assigning ranks, it calls `update_lvl_info` for each level.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.

- **`mg_load_balance_parents(mg)`:**
  - **Description:** This routine is designed for scenarios where the ranks of leaf boxes have already been determined by some other mechanism (e.g., application-specific criteria or a previous load balancing step). It then focuses on assigning ranks to the parent (non-leaf) boxes.
    1.  It iterates from `mg%highest_lvl - 1` down to `single_cpu_lvl + 1`.
    2.  For each parent box on a given level, it assigns the rank that owns the majority of its children, using `most_popular` to break ties based on the existing workload (which includes the pre-assigned leaf ranks and already processed parent ranks).
    3.  Boxes on levels at or below `single_cpu_lvl` are assigned to a `coarse_rank` determined by the most prevalent rank on `single_cpu_lvl + 1`.
    After assigning ranks, it calls `update_lvl_info` for each level.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.

### Private Helper Subroutines

- **`most_popular(list, work, n_cpu) pure integer function`:**
  - **Description:** Identifies the most frequently occurring MPI rank in the input integer `list(:)`. If there's a tie in frequency, it selects the rank among the tied ones that currently has the minimum `work` assigned (from the `work(0:n_cpu-1)` array).
  - **Arguments:** `list (integer, intent(in), :)`, `work (integer, intent(in), array)`, `n_cpu (integer, intent(in))`.

- **`update_lvl_info(mg, lvl)`:**
  - **Description:** After the `mg%boxes(:)%rank` assignments are finalized by one of the load balancing routines, this subroutine updates the per-process lists within the `mg%lvls(lvl)` structure. Specifically, for the current MPI process (identified by `mg%my_rank`), it populates:
    - `lvl%my_ids`: List of all box IDs on this `lvl` owned by the current process.
    - `lvl%my_leaves`: List of leaf box IDs on this `lvl` owned by the current process.
    - `lvl%my_parents`: List of parent box IDs on this `lvl` owned by the current process.
    - `lvl%my_ref_bnds`: List of refinement boundary box IDs on this `lvl` owned by the current process.
  - **Arguments:** `mg (type(mg_t), intent(inout))`, `lvl (type(mg_lvl_t), intent(inout))`.

## Important Variables/Constants

- **`mg%boxes(id)%rank (integer)`:** The primary output of this module; stores the MPI rank assigned to own box `id`.
- **`single_cpu_lvl (integer, local variable)`:** Calculated as `max(mg%first_normal_lvl-1, mg%lowest_lvl)`. Boxes on levels at or below `single_cpu_lvl` are typically assigned to a single MPI rank (often rank 0) because they might have varying box sizes not suitable for standard distributed communication patterns or represent too little computational work to be worth distributing.

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Absolutely fundamental. This module extensively uses `type(mg_t)`, `type(mg_box_t)`, and `type(mg_lvl_t)` and their members, including `rank`, `ids`, `leaves`, `parents`, `children`, `my_ids`, `first_normal_lvl`, `lowest_lvl`, `highest_lvl`, `n_cpu`, and `my_rank`.
- **External Libraries:**
  - None are directly invoked for MPI communication within this module's subroutines, but the logic is entirely built around an MPI-parallel environment defined by `mg%n_cpu` and `mg%my_rank`.
- **Interactions with Other Components:**
  - **`m_build_tree`:** Load balancing routines are typically called after the grid tree has been initially constructed by `m_build_tree` (which usually assigns all boxes to rank 0 by default) or after any significant structural changes to the tree (e.g., due to adaptive mesh refinement, though AMR is not explicitly detailed in `m_build_tree`).
  - **`m_allocate_storage`:** This module is called *after* load balancing. `m_allocate_storage` uses the `mg%boxes(id)%rank` assignments to allocate memory for the actual cell-centered data arrays (`mg%boxes(id)%cc`) only for those boxes assigned to the current MPI process (`mg%my_rank`).
  - **Communication Modules (e.g., `m_ghost_cells`, `m_restrict`, `m_prolong`, `m_communication`):** All inter-process communication routines rely on the `mg%boxes(id)%rank` field to determine if data needs to be sent to or received from another process, or if a neighboring box is local.
  - **Adaptive Mesh Refinement (AMR) Systems:** If the application employs AMR, strategies like `mg_load_balance_parents` become particularly important to re-balance the parent nodes after the distribution of leaf nodes changes due to refinement or coarsening.

## Usage Examples

```fortran
! Conceptual Example: Performing load balancing after tree construction.

use m_data_structures
use m_build_tree     ! For tree creation (assumed to be called prior)
use m_load_balance
use m_allocate_storage ! Typically called after load balancing

implicit none

type(mg_t) :: my_multigrid_config

! 1. Initialize the multigrid configuration and build the tree structure.
!    (This is a complex process involving other modules; mg_build_rectangle
!    usually assigns all boxes to rank 0 initially).
!    call initialize_and_build_tree(my_multigrid_config)

! 2. Perform load balancing. Choose one of the available strategies:

!    Strategy A: Simple load balancing (often suitable for uniform grids)
!    call mg_load_balance_simple(my_multigrid_config)

!    Strategy B: More advanced load balancing (tries to keep children with parents)
call mg_load_balance(my_multigrid_config)

!    Strategy C: Balance parents assuming leaves are already balanced
!    (Requires leaf ranks to be set prior to this call by another mechanism)
!    call mg_load_balance_parents(my_multigrid_config)

! 3. After load balancing, the mg%boxes(:)%rank fields are set, and the
!    mg%lvls(:)%my_ids, my_leaves, my_parents, my_ref_bnds arrays are populated.
!    Now, allocate storage for the data arrays of locally owned boxes.
call mg_allocate_storage(my_multigrid_config)

! The application can now proceed with computations, with each MPI process
! working on its assigned subset of boxes.
if (my_multigrid_config%my_rank == 0) then
    print *, "Load balancing and storage allocation complete."
end if
```
