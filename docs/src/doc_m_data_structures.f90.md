# `m_data_structures.f90`

## Overview

The `m_data_structures` module is arguably the most critical and foundational module within the `m_octree_mg` library. It serves as the central repository for all primary data structure definitions (Fortran derived types) and named constants. These definitions provide the common language and data representation used by virtually all other modules in the library to describe the multigrid hierarchy, individual grid levels, computational boxes (cells/nodes), communication buffers, boundary conditions, and various control parameters. Understanding this module is key to understanding how data is organized and manipulated throughout the octree-mg framework.

## Key Components

### Modules

- **`m_data_structures`:** Defines derived types, named constants, abstract interfaces for procedure pointers, and some utility functions related to these data structures.

### Derived Type Definitions

- **`mg_lvl_t`**: Represents data associated with a specific refinement level in the multigrid hierarchy.
  - `leaves (integer, allocatable, :)`: IDs of boxes on this level that are leaves (have no children).
  - `parents (integer, allocatable, :)`: IDs of boxes on this level that are parents (have children).
  - `ref_bnds (integer, allocatable, :)`: IDs of leaf boxes on this level that are adjacent to coarser boxes (refinement boundaries).
  - `ids (integer, allocatable, :)`: All box IDs belonging to this level across all processes.
  - `my_leaves (integer, allocatable, :)`: IDs of leaf boxes on this level owned by the current MPI process.
  - `my_parents (integer, allocatable, :)`: IDs of parent boxes on this level owned by the current MPI process.
  - `my_ref_bnds (integer, allocatable, :)`: IDs of refinement boundary boxes on this level owned by the current MPI process.
  - `my_ids (integer, allocatable, :)`: All box IDs on this level owned by the current MPI process.

- **`mg_box_t`**: Represents a single computational box (a node in the octree/quadtree).
  - `rank (integer)`: The MPI rank of the process that owns this box.
  - `id (integer)`: Unique identifier for this box (typically its index in the global `mg%boxes(:)` array).
  - `lvl (integer)`: The refinement level this box belongs to.
  - `ix(NDIM) (integer)`: Logical spatial index of this box within its level (e.g., `[i,j,k]`).
  - `parent (integer)`: ID of the parent box (or `mg_no_box`).
  - `children(2**NDIM) (integer)`: IDs of child boxes (or `mg_no_box`). `NDIM` is resolved by preprocessor.
  - `neighbors(2*NDIM) (integer)`: IDs of face-neighboring boxes (or `mg_no_box`, `mg_physical_boundary`). `NDIM` is resolved by preprocessor.
  - `r_min(NDIM) (real(dp))`: Physical coordinates of the minimum corner of this box.
  - `dr(NDIM) (real(dp))`: Grid spacing (cell size) within this box.
  - `cc (real(dp), allocatable, DTIMES(:), :)`: Cell-centered data array. Dimensions depend on `NDIM` (e.g., `(0:nc+1, 0:nc+1, 0:nc+1, mg_num_vars + n_extra_vars)` for 3D including ghost cells). The `DTIMES(:)` macro likely handles the spatial dimensions. The second dimension indexes different variables (solution, RHS, etc.).

- **`mg_buf_t`**: Represents communication buffers used for MPI data exchange, typically one per communicating pair of processes.
  - `i_send (integer)`: Current number of elements in the `send` buffer.
  - `i_recv (integer)`: Current number of elements in the `recv` buffer.
  - `i_ix (integer)`: Current number of elements in the `ix` (index/sorting key) array.
  - `ix (integer, allocatable, :)`: Array of indices or keys used for sorting data in the `send` buffer before transmission.
  - `send (real(dp), allocatable, :)`: Buffer for outgoing data.
  - `recv (real(dp), allocatable, :)`: Buffer for incoming data.

- **`mg_comm_t`**: Used to store information about communication patterns for specific operations (restriction, prolongation, ghost cells).
  - `n_send (integer, allocatable, :, :)`: Number of items to send to other processes.
  - `n_recv (integer, allocatable, :, :)`: Number of items to receive from other processes.

- **`mg_bc_t`**: Defines a boundary condition type and associated data.
  - `bc_type (integer)`: Type of boundary condition (e.g., `mg_bc_dirichlet`, `mg_bc_neumann`). Default: `mg_bc_dirichlet`.
  - `bc_value (real(dp))`: Value associated with the boundary condition (e.g., Dirichlet value, Neumann flux). Default: `0.0_dp`.
  - `boundary_cond (procedure(mg_subr_bc), pointer, nopass)`: Procedure pointer to a user-defined subroutine for applying physical boundary conditions. Default: `null()`.
  - `refinement_bnd (procedure(mg_subr_rb), pointer, nopass)`: Procedure pointer to a user-defined subroutine for handling refinement boundaries (interpolating from coarse to fine). Default: `null()`.

- **`mg_timer_t`**: Structure for simple performance timers.
  - `name (character(len=20))`: Name of the timer.
  - `t (real(dp))`: Accumulated time. Default: `0.0_dp`.
  - `t0 (real(dp))`: Start time of the current timing interval.

- **`mg_t`**: The main derived type encapsulating the entire multigrid hierarchy and its associated settings and data.
  - `tree_created (logical)`: True if the basic tree structure has been built. Default: `.false.`.
  - `is_allocated (logical)`: True if memory for box data (`cc` arrays, buffers) has been allocated. Default: `.false.`.
  - `n_extra_vars (integer)`: Number of additional user-defined cell-centered variables. Default: `0`.
  - `comm (integer)`: MPI communicator. Default: `-1`.
  - `n_cpu (integer)`: Total number of MPI processes. Default: `-1`.
  - `my_rank (integer)`: Rank of the current MPI process. Default: `-1`.
  - `box_size (integer)`: Number of cells per side for boxes at the reference level (level 1). Default: `-1`.
  - `highest_lvl (integer)`: Highest (finest) refinement level present in the grid. Default: `-1`.
  - `lowest_lvl (integer)`: Lowest (coarsest) refinement level present in the grid. Default: `-1`.
  - `first_normal_lvl (integer)`: Coarsest level at which children fully populate $2^{NDIM}$ sub-regions of their parent. Below this, parents might have only one child (special coarsening). Default: `-1`.
  - `n_boxes (integer)`: Total number of boxes in the grid across all processes. Default: `0`.
  - `box_size_lvl(mg_lvl_lo:mg_lvl_hi) (integer)`: Box size (cells per side) for each level.
  - `domain_size_lvl(NDIM, mg_lvl_lo:mg_lvl_hi) (integer)`: Effective domain size in cells for each level if uniformly refined.
  - `dr(NDIM, mg_lvl_lo:mg_lvl_hi) (real(dp))`: Grid spacing for each level.
  - `r_min(NDIM) (real(dp))`: Physical coordinates of the minimum corner of the overall domain.
  - `lvls(mg_lvl_lo:mg_lvl_hi) (type(mg_lvl_t))`: Array of level structures.
  - `boxes (type(mg_box_t), allocatable, :)`: Array containing all box structures in the grid.
  - `buf (type(mg_buf_t), allocatable, :)`: Array of communication buffers, one for each MPI process.
  - `comm_restrict (type(mg_comm_t))`: Communication pattern info for restriction.
  - `comm_prolong (type(mg_comm_t))`: Communication pattern info for prolongation.
  - `comm_ghostcell (type(mg_comm_t))`: Communication pattern info for ghost cell filling.
  - `phi_bc_data_stored (logical)`: Flag indicating if boundary condition data for the solution variable has been stored. Default: `.false.`.
  - `periodic(NDIM) (logical)`: Flags indicating periodicity in each dimension. Default: `.false.`.
  - `bc(mg_num_neighbors, mg_max_num_vars) (type(mg_bc_t))`: Array to store pre-defined boundary conditions for each face and variable.
  - `operator_type (integer)`: Type of PDE operator (e.g., `mg_laplacian`). Default: `mg_laplacian`.
  - `geometry_type (integer)`: Type of coordinate system (e.g., `mg_cartesian`). Default: `mg_cartesian`.
  - `subtract_mean (logical)`: Whether to subtract the mean from the solution (e.g., for pure Neumann Poisson). Default: `.false.`.
  - `smoother_type (integer)`: Type of multigrid smoother (e.g., `mg_smoother_gs`). Default: `mg_smoother_gs`.
  - `n_smoother_substeps (integer)`: Number of substeps for the smoother (e.g., 2 for Red-Black Gauss-Seidel). Default: `1`.
  - `n_cycle_down (integer)`: Number of smoothing steps during the downward (restriction) leg of a V-cycle. Default: `2`.
  - `n_cycle_up (integer)`: Number of smoothing steps during the upward (prolongation) leg of a V-cycle. Default: `2`.
  - `max_coarse_cycles (integer)`: Maximum number of iterations on the coarsest grid. Default: `1000`.
  - `coarsest_grid(NDIM) (integer)`: Minimum box size on the coarsest grid. Default: `2` cells per side.
  - `residual_coarse_abs (real(dp))`: Absolute residual tolerance for stopping coarse grid solve. Default: `1e-8_dp`.
  - `residual_coarse_rel (real(dp))`: Relative residual reduction factor for stopping coarse grid solve. Default: `1e-8_dp`.
  - `box_op (procedure(mg_box_op), pointer, nopass)`: Procedure pointer for the main PDE operator. Default: `null()`.
  - `box_smoother (procedure(mg_box_gsrb), pointer, nopass)`: Procedure pointer for the smoother. Default: `null()`.
  - `box_prolong (procedure(mg_box_prolong), pointer, nopass)`: Procedure pointer for the prolongation operator. Default: `null()`.
  - `n_timers (integer)`: Number of active timers. Default: `0`.
  - `timers(mg_max_timers) (type(mg_timer_t))`: Array of timer structures.

### Abstract Interfaces

These define the signatures for procedure pointers used in `mg_t` and `mg_bc_t`.
- **`mg_subr_bc`**: For user-defined physical boundary condition routines.
  - **Arguments:** `box (mg_box_t, in)`, `nc (integer, in)`, `iv (integer, in)`, `nb (integer, in)`, `bc_type (integer, out)`, `bc (real(dp), out, array)`
- **`mg_subr_rb`**: For user-defined refinement boundary routines (interpolation from coarse grid data `cgc`).
  - **Arguments:** `box (mg_box_t, inout)`, `nc (integer, in)`, `iv (integer, in)`, `nb (integer, in)`, `cgc (real(dp), in, array)`
- **`mg_box_op`**: For routines that apply the discretized PDE operator (e.g., $L\phi$).
  - **Arguments:** `mg (mg_t, inout)`, `id (integer, in)`, `nc (integer, in)`, `i_out (integer, in)`
- **`mg_box_gsrb`**: For smoother routines (e.g., Gauss-Seidel Red-Black).
  - **Arguments:** `mg (mg_t, inout)`, `id (integer, in)`, `nc (integer, in)`, `redblack_cntr (integer, in)`
- **`mg_box_prolong`**: For prolongation routines (interpolating data from a parent box to a child).
  - **Arguments:** `mg (mg_t, inout)`, `p_id (integer, in)`, `dix(NDIM) (integer, in)`, `nc (integer, in)`, `iv (integer, in)`, `fine (real(dp), out, array)`

### Important Variables/Constants (Parameters)

- **Precision:**
  - `dp (integer, parameter)`: Kind parameter for double precision reals (from `kind(0.0d0)`).
  - `i8 (integer, parameter)`: Kind parameter for 64-bit integers (from `selected_int_kind(18)`).
- **Operator Types:** `mg_laplacian`, `mg_vlaplacian` (variable coefficient), `mg_helmholtz`, `mg_vhelmholtz`, `mg_ahelmholtz` (anisotropic Helmholtz).
- **Geometry Types:** `mg_cartesian`, `mg_cylindrical`, `mg_spherical`.
- **Smoother Types:** `mg_smoother_gs` (Gauss-Seidel), `mg_smoother_gsrb` (Gauss-Seidel Red-Black), `mg_smoother_jacobi`.
- **Dimensionality:**
  - `mg_ndim (integer, parameter)`: Problem dimension, set by preprocessor `NDIM`.
- **Variable Indices (for `mg_box_t%cc` array):**
  - `mg_num_vars (integer, parameter)`: Number of predefined variables (4: phi, rhs, old_phi, residual).
  - `mg_max_num_vars (integer, parameter)`: Maximum allowed variables (10).
  - `mg_iphi (integer, parameter)`: Index for the solution variable $\phi$.
  - `mg_irhs (integer, parameter)`: Index for the right-hand side $f$.
  - `mg_iold (integer, parameter)`: Index for the previous solution (used in correction scheme).
  - `mg_ires (integer, parameter)`: Index for the residual $f - L\phi$.
  - `mg_iveps (integer, parameter)`: Index for isotropic variable coefficient $\epsilon$.
  - `mg_iveps1`, `mg_iveps2`, `mg_iveps3 (integer, parameter)`: Indices for anisotropic coefficients $\epsilon_x, \epsilon_y, \epsilon_z$ (existence depends on `NDIM`).
- **Grid Level Limits:** `mg_lvl_lo`, `mg_lvl_hi` (parameters defining min/max level indices).
- **Boundary Condition Types:** `mg_bc_dirichlet`, `mg_bc_neumann`, `mg_bc_continuous`.
- **Special Box/Neighbor IDs:**
  - `mg_no_box (integer, parameter)`: Indicates no box (e.g., a child pointer for a leaf box).
  - `mg_physical_boundary (integer, parameter)`: Indicates a neighbor link pointing to a physical domain boundary.
- **Timers:** `mg_max_timers (integer, parameter)`: Maximum number of timers.
- **Geometric/Topological Constants (NDIM-dependent, defined via C-preprocessor):**
  - `mg_num_children (integer, parameter)`: $2^{NDIM}$.
  - `mg_child_dix(NDIM, mg_num_children) (integer, parameter)`: Index offsets of children relative to parent's logical index space.
  - `mg_child_rev(mg_num_children, NDIM) (integer, parameter)`: Mapping to find a child's index if looking from an adjacent parent.
  - `mg_child_adj_nb(mg_num_children/2, 2*NDIM) (integer, parameter)`: Children adjacent to a specific neighbor face.
  - `mg_child_low(NDIM, mg_num_children) (logical, parameter)`: Flags indicating if a child is on the "low" side in each dimension.
  - `mg_num_neighbors (integer, parameter)`: $2 \times NDIM$.
  - `mg_neighb_lowx`, `mg_neighb_highx`, etc. (integer, parameter): Named constants for neighbor directions.
  - `mg_neighb_dix(NDIM, mg_num_neighbors) (integer, parameter)`: Index offsets to find neighbors.
  - `mg_neighb_low(mg_num_neighbors) (logical, parameter)`: Flags indicating if a neighbor is in the "low" direction.
  - `mg_neighb_high_pm(mg_num_neighbors) (integer, parameter)`: `[-1, 1]` representation for low/high neighbors.
  - `mg_neighb_rev(mg_num_neighbors) (integer, parameter)`: Opposite direction for each neighbor.
  - `mg_neighb_dim(mg_num_neighbors) (integer, parameter)`: Dimension corresponding to each neighbor direction.

### Module Subroutines/Functions

Besides type/constant definitions, the module also contains several public utility functions:
- **`mg_has_children(box)`**: Elemental function, returns `.true.` if `box%children(1)` is not `mg_no_box`.
- **`mg_ix_to_ichild(ix)`**: Calculates a child's index (1 to $2^{NDIM}$) within its parent's `children` array based on its logical spatial index `ix`.
- **`mg_get_child_offset(mg, id)`**: Returns the spatial offset of a child box relative to its parent's origin, in units of child cell width.
- **`mg_highest_uniform_lvl(mg)`**: Determines the highest level up to which the grid is uniformly refined.
- **`mg_number_of_unknowns(mg)`**: Calculates the total number of unknowns (degrees of freedom) on all leaf boxes.
- **`mg_get_face_coords(box, nb, nc, x)`**: Calculates physical coordinates of cell centers on a specific face `nb` of a `box`.
- **Timer routines:** `mg_add_timer`, `mg_timer_start`, `mg_timer_end`, `mg_timers_show` for managing and displaying simple performance timers.

## Usage Examples

```fortran
! Example: Declaring mg_t and accessing box data
use m_data_structures
implicit none

type(mg_t) :: my_multigrid_solver
type(mg_box_t) :: current_box
integer :: box_id, level_num, num_cells_per_side
real(dp) :: solution_value, rhs_value
integer :: i, j, k ! Cell indices

! ... (Assume my_multigrid_solver is initialized, tree built, and storage allocated)
! ... (Assume box_id is a valid ID for a box owned by the current process)

! Get a reference to a specific box (conceptual, direct access might be via mg%boxes(box_id))
! current_box = my_multigrid_solver%boxes(box_id) ! This is how you'd actually get it

! For demonstration, let's assume current_box is somehow populated:
! (This is just for showing access, not functional code for getting a box)
level_num = 1
current_box%lvl = level_num
current_box%id = box_id
num_cells_per_side = my_multigrid_solver%box_size_lvl(level_num)

! Allocate cell-centered data for the box (nc = num_cells_per_side)
! The DTIMES macro handles dimensions, assuming 0:nc+1 for ghost cells
! Let's say mg_num_vars + n_extra_vars = 5
#if NDIM == 3
  allocate(current_box%cc(0:num_cells_per_side+1, 0:num_cells_per_side+1, 0:num_cells_per_side+1, 5))
#elif NDIM == 2
  allocate(current_box%cc(0:num_cells_per_side+1, 0:num_cells_per_side+1, 5))
#else
  allocate(current_box%cc(0:num_cells_per_side+1, 5))
#endif
current_box%cc = 0.0_dp ! Initialize

! Example: Accessing data within a 3D box at cell (i,j,k)
i = num_cells_per_side / 2
j = num_cells_per_side / 2
k = num_cells_per_side / 2

! Set and get solution (phi) and RHS values
#if NDIM == 3
  current_box%cc(i,j,k, mg_iphi) = 1.23_dp
  current_box%cc(i,j,k, mg_irhs) = 4.56_dp
  solution_value = current_box%cc(i,j,k, mg_iphi)
  rhs_value      = current_box%cc(i,j,k, mg_irhs)
  print *, "3D Example: Solution at box", box_id, "cell (",i,j,k,") is", solution_value
#elif NDIM == 2
  current_box%cc(i,j, mg_iphi) = 1.23_dp
  current_box%cc(i,j, mg_irhs) = 4.56_dp
  solution_value = current_box%cc(i,j, mg_iphi)
  rhs_value      = current_box%cc(i,j, mg_irhs)
  print *, "2D Example: Solution at box", box_id, "cell (",i,j,") is", solution_value
#else
  current_box%cc(i, mg_iphi) = 1.23_dp
  current_box%cc(i, mg_irhs) = 4.56_dp
  solution_value = current_box%cc(i, mg_iphi)
  rhs_value      = current_box%cc(i, mg_irhs)
  print *, "1D Example: Solution at box", box_id, "cell (",i,") is", solution_value
#endif

! Accessing other properties of the mg_t structure
print *, "Multigrid solver configured for operator type:", my_multigrid_solver%operator_type
print *, "Grid periodicity in X-dim:", my_multigrid_solver%periodic(1)

if (allocated(current_box%cc)) deallocate(current_box%cc)

```

## Dependencies and Interactions

- **Internal Dependencies:**
  - Primarily uses intrinsic Fortran capabilities.
  - Uses `mpi` module in timer functions (`mg_timer_start`, `mg_timer_end`, `mg_timers_show`) for `mpi_wtime()` and `mpi_reduce()`.
- **External Libraries:**
  - `cpp_macros.h`: This C-preprocessor include is used for `NDIM` (which defines the dimensionality compiled into the library) and `DTIMES` (which likely expands to the appropriate number of colons for array slicing based on `NDIM`).
- **Interactions with Other Components:**
  - This module is fundamental to the entire `m_octree_mg` library. Nearly every other Fortran module in the `src/` directory will `use m_data_structures` to access the definitions of `mg_t`, `mg_box_t`, `mg_lvl_t`, various constants (like `mg_iphi`, `mg_bc_dirichlet`), and procedure interfaces. It provides the common vocabulary and data framework for all other components that build, manipulate, solve on, or communicate about the multigrid structure.
```
