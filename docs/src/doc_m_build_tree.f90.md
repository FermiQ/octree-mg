# `m_build_tree.f90`

## Overview

The `m_build_tree` module is a cornerstone of the `m_octree_mg` library, tasked with the construction and structural definition of the hierarchical grid. This grid is typically a quadtree in two dimensions or an octree in three dimensions. The module handles the logic for creating the initial coarse grid based on specified domain parameters and then recursively building finer grid levels. This process involves instantiating box structures (nodes in the tree), defining their properties (level, spatial index, physical coordinates), establishing parent-child relationships between boxes at different levels, and setting up initial neighbor connectivity.

## Key Components

### Modules

- **`m_build_tree`:** This module encapsulates all subroutines necessary for the initial generation and hierarchical organization of the multigrid tree structure.

### Functions/Subroutines

- **`mg_build_rectangle(mg, domain_size, box_size, dx, r_min, periodic, n_finer)`:**
  - **Description:** This is the main public subroutine for constructing the complete multi-level grid hierarchy. It starts by defining the properties of the base level (level 1, typically the finest regular grid) such as grid spacing (`dr`), box size, and domain size. It then determines the parameters for coarser levels by either doubling the grid spacing or halving the number of boxes per dimension, until a specified coarsest grid criterion (`mg%coarsest_grid`) is met. It also defines metadata for finer levels if applicable. The subroutine allocates the main `mg%boxes` array to hold all box structures. It then explicitly creates all boxes for the determined `mg%lowest_lvl`, setting their geometric properties, initial neighbor information (handling periodic and physical boundaries), and parent/child links (initially `mg_no_box`). Subsequently, it iteratively calls `mg_add_children` (or `add_single_child` for specific coarsening strategies) to populate finer levels and uses helper subroutines (`mg_set_leaves_parents`, `mg_set_next_level_ids`, `mg_set_neighbors_lvl`) to complete the tree structure by identifying parent/leaf nodes at each level and ensuring full neighbor connectivity. Finally, it marks the tree as created (`mg%tree_created = .true.`).
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The main multigrid data structure to be populated.
    - `domain_size (integer(NDIM), intent(in))`: Total number of cells in each dimension for the base level.
    - `box_size (integer, intent(in))`: Number of cells per side for each box at the base level.
    - `dx (real(dp)(NDIM), intent(in))`: Grid spacing in each dimension at the base level.
    - `r_min (real(dp)(NDIM), intent(in))`: Physical coordinate of the domain's minimum corner.
    - `periodic (logical(NDIM), intent(in))`: Flags indicating periodicity in each dimension.
    - `n_finer (integer, intent(in))`: Number of additional box structures to allocate space for, anticipating future refinements.

- **`mg_add_children(mg, id)`:**
  - **Description:** Creates a full set of child boxes (e.g., 4 for `NDIM=2`, 8 for `NDIM=3`) for a given parent box `id`. It increments `mg%n_boxes`, assigns IDs to the new children, and links them to the parent `id`. Properties for each child (rank, logical index `ix`, level, parent pointer, initial neighbors, physical minimum `r_min`, and grid spacing `dr`) are derived from the parent box and the grid hierarchy. Boundary conditions for children facing the domain edge are set based on the parent's boundary conditions.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `id (integer, intent(in))`: The ID of the parent box to which children will be added.

- **`mg_set_leaves_parents(boxes, level)`:**
  - **Description:** Scans all boxes listed in `level%ids`. For each box, it checks if it has children. If it does, the box ID is added to `level%parents`. If it does not, the ID is added to `level%leaves`. This populates lists of parent nodes and terminal (leaf) nodes at a specific grid level.
  - **Arguments:**
    - `boxes (type(mg_box_t), intent(in))`: The array of all box structures.
    - `level (type(mg_lvl_t), intent(inout))`: The level structure whose `parents` and `leaves` lists are to be populated.

- **`mg_set_next_level_ids(mg, lvl)`:**
  - **Description:** Populates the `mg%lvls(lvl+1)%ids` array. It collects all the children from the parent boxes found in `mg%lvls(lvl)%parents` and lists their IDs, effectively defining all boxes that exist at level `lvl+1`.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `lvl (integer, intent(in))`: The current level whose parents' children will form the next level.

- **`mg_set_refinement_boundaries(boxes, level)`:**
  - **Description:** Identifies leaf boxes on a given `level` that are adjacent to coarser boxes (i.e., a neighbor of the leaf box is a parent box on the same level). Such leaf boxes form the finer side of a refinement boundary. Their IDs are stored in `level%ref_bnds`.
  - **Arguments:**
    - `boxes (type(mg_box_t), intent(in))`: The array of all box structures.
    - `level (type(mg_lvl_t), intent(inout))`: The level structure for which refinement boundaries are identified.

- **`mg_set_neighbors_lvl(mg, lvl)`:**
  - **Description:** Ensures all boxes on a specified `lvl` have their neighbor information correctly and fully established. It iterates through `mg%lvls(lvl)%ids` and calls the private `set_neighbs` routine for each box.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure.
    - `lvl (integer, intent(in))`: The level for which neighbor information is to be finalized.

### Private Helper Subroutines (Notable)
- `add_single_child(mg, id, n_boxes_lvl)`: Creates a single child for a box, used in specific scenarios where the effective grid resolution changes between levels differently than the standard refinement.
- `set_neighbs(boxes, id)`: Orchestrates neighbor finding for a single box `id` by calling `find_neighb` for each neighbor direction if not already set.
- `find_neighb(boxes, id, nb)`: A key routine that determines the ID of a neighbor for box `id` in direction `nb`. It might navigate up to the parent level and across to the parent's neighbor, then down to the appropriate child on the adjacent branch of the tree.

## Important Variables/Constants

The tree construction process is heavily influenced by:
- Inputs to `mg_build_rectangle`: `domain_size`, `box_size`, `dx`, `r_min`, `periodic`, `n_finer`.
- Parameters within the `mg_t` structure (often set by `mg_init_default` or other setup routines):
  - `mg%coarsest_grid`: Defines a limit on how coarse the grid hierarchy can be.
  - `mg%smoother_type`: Can influence coarsening strategy, especially `mg_smoother_gs`.
- Constants from `m_data_structures`:
  - `NDIM`: Number of spatial dimensions.
  - `mg_lvl_lo`, `mg_lvl_hi`: Bounds for multigrid level indexing.
  - `mg_no_box`, `mg_physical_boundary`: Special values for box IDs/neighbor links.
  - `mg_num_children`, `mg_num_neighbors`: Geometric constants (2/4 in 1D, 4/8 in 2D, 8/26 in 3D, though neighbors are typically face-neighbors).
  - Various `mg_child_...` arrays/functions: Provide geometric relationships between parent and child boxes.

## Usage Examples

```fortran
! Conceptual example of building the initial tree:
use m_data_structures
use m_build_tree
! Potentially use m_default_settings to initialize mg_t
! use m_default_settings, only: mg_init_default

type(mg_t) :: my_grid
integer    :: domain_cells(2), cells_per_box_side
real(dp)   :: cell_spacing(2), domain_origin(2)
logical    :: domain_is_periodic(2)
integer    :: space_for_future_boxes

! call mg_init_default(my_grid) ! Initialize my_grid with library defaults

! Define parameters for a 2D grid
domain_cells         = [256, 256]  ! Total cells in underlying finest grid
cells_per_box_side   = 32          ! Cells per side for each box
cell_spacing         = [1.0_dp/256.0_dp, 1.0_dp/256.0_dp]
domain_origin        = [0.0_dp, 0.0_dp]
domain_is_periodic   = [.false., .false.]
space_for_future_boxes = 2000 ! Buffer for boxes added by dynamic refinement

! Construct the tree
call mg_build_rectangle(my_grid, domain_cells, cells_per_box_side, &
                         cell_spacing, domain_origin, domain_is_periodic, &
                         space_for_future_boxes)

if (my_grid%tree_created) then
  print *, "Octree/Quadtree structure built. Total boxes created:", my_grid%n_boxes
  print *, "Lowest level:", my_grid%lowest_lvl, "Highest level:", my_grid%highest_lvl
  ! Next steps typically include:
  ! - Distributing boxes among MPI processes (m_load_balance)
  ! - Allocating memory for cell data (m_allocate_storage)
else
  print *, "Failed to build the tree structure."
end if

! Note: Direct public subroutines for dynamic refinement (refine_box) or
! coarsening (coarsen_box) are not part of this specific module.
! However, mg_add_children is a building block that would be used by such features.
```

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** This is the most critical dependency. The module extensively uses and manipulates the `mg_t`, `mg_box_t`, and `mg_lvl_t` derived types and their numerous members (e.g., `boxes`, `lvls`, `n_boxes`, `tree_created`, `ix`, `lvl`, `parent`, `children`, `neighbors`, `r_min`, `dr`, `ids`, `leaves`, `parents`, `ref_bnds`). It also relies on various constants defined in `m_data_structures` that relate to grid geometry and special marker values.
- **External Libraries:**
  - **`cpp_macros.h`:** Used for C-preprocessor macros such as `KJI_DO_VEC` and `CLOSE_DO`, which are likely shorthand for multi-dimensional loop constructs.
- **Interactions with Other Components:**
  - **`m_allocate_storage`:** After `mg_build_rectangle` successfully executes (setting `mg%tree_created = .true.`), the `m_allocate_storage` module is typically called to allocate memory for the actual cell-centered data (`cc` arrays) within each box (`mg%boxes(id)%cc`) and for communication buffers. The `mg_build_rectangle` subroutine itself allocates the `mg%boxes` array (array of structures).
  - **`m_load_balance`:** The initial tree built by `mg_build_rectangle` assigns all boxes to rank 0. A load balancing module (`m_load_balance`) is essential in a parallel environment to distribute these boxes among available MPI processes. Dynamic tree modifications (refinement/coarsening) would also likely trigger load balancing.
  - **Initialization Routines (e.g., `m_default_settings`):** Modules that set default values or overall parameters for the `mg_t` type are usually called before `mg_build_rectangle`.
  - **Boundary Condition Modules:** The `periodic` flags passed to `mg_build_rectangle` and the assignment of `mg_physical_boundary` to neighbor links directly influence how boundary conditions are identified and applied by other parts of the solver.
  - **Refinement/Coarsening Modules:** While not fully exposed as public interfaces in this file, any higher-level modules that implement adaptive mesh refinement (AMR) would use routines like `mg_add_children` (for refinement) and would need mechanisms to remove or deactivate boxes (for coarsening), which would modify the tree structure created by this module.
```
