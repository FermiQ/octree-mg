# Overview of the `single_module` Directory

## Overview

The `single_module` directory within the `octree-mg` library distribution offers a convenient way to integrate the core functionalities of the library into other Fortran projects. Instead of managing multiple source files and module dependencies from the main `src/` directory, this directory provides options to include the entire library (with some considerations for the free-space Poisson solver) as a single Fortran module file. This approach can significantly simplify the build process for projects that utilize `octree-mg`.

The primary mechanism involves a main source file, `m_octree_mg.f90`, which uses a C-preprocessor directive (`NDIM`) to generate dimension-specific code (1D, 2D, or 3D). Pre-generated versions for each dimension are also typically provided.

## Key Files and Their Purpose

-   **`m_octree_mg.f90`**:
    -   This is the main, preprocessor-driven source file. It is designed to be compiled with a flag like `-DNDIM=1`, `-DNDIM=2`, or `-DNDIM=3`.
    -   When compiled, it effectively concatenates all the core modules from the `src/` directory (e.g., `m_data_structures`, `m_build_tree`, `m_load_balance`, `m_ghost_cells`, `m_allocate_storage`, `m_restrict`, `m_communication`, `m_prolong`, `m_multigrid`, and the various operator modules like `m_laplacian`, `m_helmholtz`, etc.) into a single dimension-specific module.
    -   The resulting Fortran module name will depend on the `NDIM` definition (e.g., `m_octree_mg_1d`, `m_octree_mg_2d`, `m_octree_mg_3d`).

-   **`m_octree_mg_1d.f90`**:
    -   A pre-generated single-file Fortran module for 1D simulations.
    -   This is likely the output of compiling `single_module/m_octree_mg.f90` with `-DNDIM=1`.
    -   It defines a Fortran module named `m_octree_mg_1d`.

-   **`m_octree_mg_2d.f90`**:
    -   A pre-generated single-file Fortran module for 2D simulations.
    -   Likely the output of compiling `single_module/m_octree_mg.f90` with `-DNDIM=2`.
    -   It defines a Fortran module named `m_octree_mg_2d`.

-   **`m_octree_mg_3d.f90`**:
    -   A pre-generated single-file Fortran module for 3D simulations.
    -   Likely the output of compiling `single_module/m_octree_mg.f90` with `-DNDIM=3`.
    -   It defines a Fortran module named `m_octree_mg_3d`.

-   **`to_single_module.sh`**:
    -   A shell script responsible for generating the dimension-specific single-file modules (`m_octree_mg_1d.f90`, `m_octree_mg_2d.f90`, `m_octree_mg_3d.f90`) from the main `m_octree_mg.f90` source file and, by extension, from the individual modules in the `src/` directory. It likely invokes the Fortran compiler with the appropriate preprocessor flags.

-   **`m_free_space.f90`**:
    -   This file, also present in the `single_module` directory (and originally in `src/`), provides the `m_free_space` module.
    -   This module implements a specialized Poisson solver for 3D free-space (open) boundary conditions.
    -   It is intended to be used alongside `m_octree_mg_3d.f90` if this specific boundary condition handling is required. It is separate because its primary dependency, the `Poisson_Solver` FFT library, is a more complex external dependency.

## How to Use

Developers wishing to integrate `octree-mg` into their projects using the single-module approach have two primary options:

1.  **Use Pre-generated Files:**
    -   Copy the appropriate dimension-specific file (`m_octree_mg_1d.f90`, `m_octree_mg_2d.f90`, or `m_octree_mg_3d.f90`) directly into their project's source directory.
    -   Compile this file as part of their project.
    -   In their own Fortran code, `use` the corresponding module (e.g., `use m_octree_mg_2d`).

2.  **Compile `m_octree_mg.f90` Manually:**
    -   Copy `single_module/m_octree_mg.f90` into their project.
    -   Compile this file using their Fortran compiler, ensuring the C preprocessor is invoked and the `NDIM` macro is defined (e.g., `gfortran -cpp -DNDIM=3 m_octree_mg.f90 -c`).
    -   The compiled module object file can then be linked, and the corresponding dimension-specific module (e.g., `m_octree_mg_3d`) can be `use`d.

By `use`ing one of these generated modules (e.g., `use m_octree_mg_3d`), all the public entities (derived types like `mg_t`, subroutines like `mg_build_rectangle`, `mg_fas_fmg`, constants like `mg_iphi`) from the original individual `src/` modules become available to the user's code.

## Note on `m_free_space.f90`

-   If the specialized 3D free-space Poisson solver is required, the `single_module/m_free_space.f90` file should also be compiled into the project.
-   The user's code would then also `use m_free_space`.
-   The `m_free_space` module itself internally `use`s `m_octree_mg_3d` (when compiled as part of the single-module system or if that module is otherwise available), indicating it's designed to work with the 3D version of the core library.

## Redundancy with `docs/src/` Documentation

The detailed documentation for the individual components, data structures, and algorithms that are bundled into these single-module files (e.g., `m_octree_mg_1d.f90`, `m_octree_mg_2d.f90`, `m_octree_mg_3d.f90`) can be found in the `docs/src/` directory. The single-module files are primarily concatenations and preprocessor adaptations of these original source files. Therefore, for understanding the specifics of `type(mg_t)`, `mg_laplacian`, `mg_multigrid`, etc., users should refer to the respective `doc_m_*.f90.md` files.

## Exclusion of `poisson_3d_fft`

It is important to note that the `poisson_3d_fft` library, which is a dependency for the `m_free_space` module, is **not** directly bundled into the `m_octree_mg_Xd.f90` single-module files. If `m_free_space` is used, the `poisson_3d_fft` library (and its own dependencies, like an MPI library and potentially ABINIT's XC routines if XC support is enabled) must be compiled and linked separately. The `single_module/m_free_space.f90` file only contains the interface code to `poisson_3d_fft`, not the solver itself.
