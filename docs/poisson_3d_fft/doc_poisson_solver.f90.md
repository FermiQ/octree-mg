# `poisson_solver.f90` (Module `Poisson_Solver`)

## Overview

The `Poisson_Solver` module is a comprehensive package designed to solve the Poisson equation, typically $\nabla^2 V(\mathbf{x}) = -4\pi \rho(\mathbf{x})$ (though the $4\pi$ factor might be handled by scaling within the solver or kernel), using Fast Fourier Transform (FFT) methods in three dimensions. This solver is versatile, supporting various boundary conditions crucial for different physical scenarios:
-   **'F' (Free):** For isolated systems where the potential decays to zero at infinity.
-   **'S' (Surface):** For systems periodic in two dimensions (xz-plane) but isolated in the third (y-direction).
-   **'P' (Periodic):** For systems periodic in all three dimensions.

The module is structured as a primary Fortran module file (`poisson_solver.f90`) that includes several other `.f90` files, each encapsulating specific parts of the solver's functionality, such as kernel generation, the main FFT solution process, exchange-correlation (XC) calculations, gradient computations, low-level FFTs, and the generation of interpolating scaling functions. It also integrates MPI for parallel execution and data distribution. This solver is notably used by the `m_free_space` module within the `m_octree_mg` library to handle open boundary conditions.

## Key Components

### Module Definition
-   **`Poisson_Solver`**: The main module a user interacts with. It re-exports public entities from the included files.

### Public Module Variable
-   **`iproc_verbose (integer)`**: A control variable for verbosity. If set to `0` (typically for MPI rank 0), additional information about the solver's execution may be printed. Default is `-1` (no verbose output).

### Public Subroutines

-   **`createKernel(geocode, n01, n02, n03, hx, hy, hz, itype_scf, iproc, nproc, kernel)`**:
    -   **Description:** This subroutine is responsible for allocating and computing the Poisson kernel (Green's function) in reciprocal (Fourier) space. The kernel's form depends on the specified boundary conditions (`geocode`) and the order of interpolating scaling functions (`itype_scf`) used for its construction.
        - For 'P' (Periodic BCs), the kernel is analytical ($1/k^2$) and this routine might perform minimal setup.
        - For 'F' (Free BCs) and 'S' (Surface BCs), the construction is more involved, often utilizing interpolating scaling functions to represent the Green's function and Gaussian quadrature for related integrals. This logic is primarily found within `build_kernel.f90`.
    -   **Arguments (Key):**
        -   `geocode (character(len=1), intent(in))`: Boundary condition type ('F', 'S', or 'P').
        -   `n01, n02, n03 (integer, intent(in))`: Global dimensions of the real-space grid.
        -   `hx, hy, hz (real(kind=8), intent(in))`: Grid spacings.
        -   `itype_scf (integer, intent(in))`: Order of interpolating scaling functions.
        -   `iproc, nproc (integer, intent(in))`: MPI rank and total number of processes.
        -   `kernel (real(kind=8), pointer, intent(out))`: Output pointer to the allocated and computed kernel array in reciprocal space.
    -   **Primary Included File:** `build_kernel.f90`.

-   **`PSolver(geocode, datacode, iproc, nproc, n01, n02, n03, ixc, hx, hy, hz, rhopot, karray, pot_ion, eh, exc, vxc, offset, sumpion, nspin)`**:
    -   **Description:** This is the main solver routine. Given a density distribution $\rho$ (in `rhopot` on input), it computes the electrostatic potential $V$ by:
        1.  Performing a forward 3D FFT on $\rho$.
        2.  Multiplying the result by the pre-computed `karray` (the Green's function $G(k)$) in reciprocal space.
        3.  Performing an inverse 3D FFT to obtain $V(\mathbf{x})$.
        The routine handles different data distributions (`datacode`: 'G' for global, 'D' for distributed data across MPI processes) and boundary conditions (`geocode`). It can optionally compute exchange-correlation (XC) energies and potentials if `ixc /= 0` (interfacing with ABINIT's XC library) and can add an external ionic potential (`pot_ion`) to the result.
    -   **Arguments (Key):**
        -   `geocode, datacode (character(len=1), intent(in))`: Boundary and data distribution codes.
        -   `iproc, nproc, n01, n02, n03 (integer, intent(in))`: MPI and grid dimension parameters.
        -   `ixc (integer, intent(in))`: Exchange-correlation functional code (0 for none).
        -   `hx, hy, hz (real(kind=8), intent(in))`: Grid spacings.
        -   `rhopot (real(kind=8), dimension(*), intent(inout))`: On input, contains the density $\rho$. On output, it's overwritten with the computed potential $V_H (+ V_{XC} + V_{ion}$ if applicable).
        -   `karray (real(kind=8), dimension(*), intent(in))`: The pre-computed Poisson kernel from `createKernel`.
        -   `pot_ion (real(kind=8), dimension(*), intent(inout))`: External ionic potential / XC potential output.
        -   `eh, exc, vxc (real(kind=8), intent(out))`: Hartree energy, XC energy, and $\int \rho V_{XC} d\mathbf{x}$.
        -   `offset (real(kind=8), intent(in))`: Potential offset for periodic calculations.
        -   `sumpion (logical, intent(in))`: If true and `ixc /=0`, adds `pot_ion` to `rhopot`.
        -   `nspin (integer, intent(in))`: Number of spin components (1 or 2).
    -   **Primary Included Files:** `psolver_main.f90` (defines `PSolver`), `psolver_base.f90` (specific P/S/F solvers), `xcenergy.f90`, `3dgradient.f90`, `fft3d.f90`.

-   **`PS_dim4allocation(geocode, datacode, iproc, nproc, n01, n02, n03, ixc, n3d, n3p, n3pi, i3xcsh, i3s)`**:
    -   **Description:** Calculates the necessary array dimensions for `rhopot` and `pot_ion` based on the problem setup (boundary conditions, data distribution, XC functional, MPI configuration). This is crucial for users to correctly allocate memory before calling `PSolver`, especially when data is distributed or XC calculations require larger ghost regions for gradients.
    -   **Arguments (Key outputs):**
        -   `n3d, n3p, n3pi (integer, intent(out))`: Calculated third dimensions for density, potential, and ionic potential arrays for the current process.
        -   `i3xcsh, i3s (integer, intent(out))`: Shifts and starting indices for array slicing, especially relevant for GGA XC calculations with distributed data.
    -   **Primary Included File:** `psolver_main.f90` (which calls `P/S/F_FFT_dimensions`).

-   **`P_FFT_dimensions`, `S_FFT_dimensions`, `F_FFT_dimensions`**:
    -   **Description:** Helper subroutines, likely made public through `PS_dim4allocation` or available directly, that determine the padded FFT grid dimensions (`n1, n2, n3`), real-space data dimensions (`md1, md2, md3`), and reciprocal-space kernel dimensions (`nd1, nd2, nd3`) required for the specified boundary conditions ('P', 'S', or 'F') and original grid sizes (`n01, n02, n03`). These ensure that array dimensions are compatible with FFT algorithms (e.g., products of small primes) and MPI data distribution.
    -   **Primary Included File:** `psolver_main.f90`.

### Role of Included Files

The `Poisson_Solver` module achieves its functionality by including several source files:
-   **`psolver_main.f90`**: Contains the main `PSolver` routine, which acts as a dispatcher to boundary-condition-specific routines, and the `PS_dim4allocation` routine with its helpers (`P/S/F_FFT_dimensions`).
-   **`build_kernel.f90`**: Contains `createKernel` and its helper subroutines (`Free_Kernel`, `Surfaces_Kernel`, `calculates_green_opt`, `calculates_green_opt_muzero`, `gequad`, `kernelfft`, `inserthalf`). It is responsible for constructing the $G(k)$ kernel, often using interpolating scaling functions (see `scaling_function.f90`) for 'F' and 'S' boundary conditions.
-   **`psolver_base.f90`**: Contains the core FFT-based solution algorithms for each boundary condition type: `P_PoissonSolver`, `S_PoissonSolver`, and `F_PoissonSolver`. These routines handle the forward FFT of density, multiplication by the kernel in Fourier space, and the backward FFT to get the potential. They also manage data transpositions and MPI communication patterns specific to the 3D FFTs.
-   **`xcenergy.f90`**: Contains `xc_energy` which interfaces with an external XC library (like ABINIT's) to compute exchange-correlation energy and potential if `ixc /= 0`. It also includes `vxcpostprocessing` for applying White-Bird corrections for GGA functionals.
-   **`3dgradient.f90`**: Contains `calc_gradient` to compute density gradients needed for GGA XC functionals, and `wb_correction` which is part of the White-Bird correction scheme.
-   **`fft3d.f90`**: Provides the low-level 1D FFT routine (`fftstp`) used by the 3D FFTs, and `ctrig` to precompute trigonometric factors. It also includes `fourier_dim` to determine FFT-friendly dimensions.
-   **`scaling_function.f90`**: Contains routines like `scaling_function`, `wavelet_function`, and various `*_trans_*` and `scf_recursion_*` subroutines. These are used to compute the values of interpolating scaling functions of different orders, which are then used in `build_kernel.f90` for constructing accurate Green's functions for non-periodic boundary conditions.

## Important Variables/Constants

-   **`iproc_verbose (integer, public)`**: Module variable to control verbose output (typically rank 0 prints if `iproc_verbose == 0`).
-   **`geocode (character(len=1))`**: Input to `PSolver` and `createKernel`, specifies boundary conditions: 'F' (Free), 'S' (Surface), 'P' (Periodic).
-   **`datacode (character(len=1))`**: Input to `PSolver`, specifies data layout: 'G' (Global - all data on all procs), 'D' (Distributed).
-   **`ixc (integer)`**: Input to `PSolver`, specifies the exchange-correlation functional to use (ABINIT convention, 0 means no XC).
-   **`itype_scf (integer)`**: Input to `createKernel`, defines the order of interpolating scaling functions used in kernel construction (e.g., 8, 14, 16).

## Dependencies and Interactions

-   **External Dependencies:**
    -   **MPI:** Heavily used throughout for parallel data distribution, FFTs, and global reductions (e.g., `MPI_ALLTOALL`, `MPI_ALLGATHERV`, `MPI_ALLREDUCE`, `MPI_BCAST`).
    -   **ABINIT Modules (referenced in comments):** `defs_xc`, `defs_basis`, `defs_datatypes`. These are required if `ixc /= 0` for exchange-correlation calculations. They are not part of the `poisson_3d_fft` package itself.
-   **Internal Dependencies (Inter-file within the module):**
    -   `poisson_solver.f90` acts as the main include file.
    -   `psolver_main.f90` calls routines in `psolver_base.f90` and `xcenergy.f90`.
    -   `build_kernel.f90` calls routines in `scaling_function.f90` and `fft3d.f90` (for `kernelfft`).
    -   `psolver_base.f90` calls routines in `fft3d.f90`.
    -   `xcenergy.f90` calls `drivexc` (external ABINIT interface) and routines in `3dgradient.f90`.
-   **Interaction with `m_octree_mg` library:**
    -   The `m_free_space` module of the `m_octree_mg` library uses `Poisson_Solver`'s `createKernel` and `PSolver` subroutines. This is typically done to solve the Poisson equation on the coarsest grid level of the multigrid hierarchy when free-space (open) boundary conditions are needed for the overall problem. The FFT solver provides the global solution component, which then informs the boundary conditions for the multigrid solver on finer levels.

## Usage Examples

The "QUICK INSTRUCTION FOR THE IMPATIENT" in `poisson_solver.f90` provides a good starting point:

```fortran
! Conceptual Example: Solving Poisson's equation using Poisson_Solver
! (Assumes serial execution for simplicity in this example, but designed for MPI)

module MyPoissonProblem
  use Poisson_Solver
  use mpi ! Poisson_Solver itself uses MPI, so an MPI environment is expected.
  implicit none

  subroutine solve_density_for_potential()
    real(kind=8), pointer :: kernel_array(:)
    real(kind=8), dimension(:,:,:), allocatable :: rhopot_data
    real(kind=8) :: hartree_energy, hx_spacing, hy_spacing, hz_spacing
    real(kind=8) :: dummy_exc, dummy_vxc, offset_val
    integer :: nx_dim, ny_dim, nz_dim, ierr
    integer :: my_rank, num_procs
    character(len=1) :: geometry_code
    logical :: sum_ionic_potential
    integer :: xc_functional_code, spin_components

    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    ! Define grid parameters
    nx_dim = 64; ny_dim = 64; nz_dim = 64  ! Example dimensions
    hx_spacing = 0.5d0; hy_spacing = 0.5d0; hz_spacing = 0.5d0
    geometry_code = 'F' ! Free boundary conditions

    ! Allocate and initialize density data (rhopot_data)
    ! For 'G' (Global) datacode, full array on all procs.
    ! For 'D' (Distributed), use PS_dim4allocation to get local sizes.
    ! For simplicity, assuming global allocation here.
    allocate(rhopot_data(nx_dim, ny_dim, nz_dim))
    ! ... Fill rhopot_data with density values ...
    rhopot_data = 1.0d0 ! Example: uniform density

    ! 1. Create the Poisson Kernel
    !    itype_scf=14 is a common choice mentioned.
    !    (kernel_array is a pointer, createKernel will allocate it)
    call createKernel(geometry_code, nx_dim, ny_dim, nz_dim, &
                      hx_spacing, hy_spacing, hz_spacing, &
                      14, my_rank, num_procs, kernel_array)

    ! 2. Solve for the potential
    xc_functional_code = 0  ! No XC calculation
    offset_val = 0.0d0      ! No offset for Free BC
    sum_ionic_potential = .false.
    spin_components = 1

    ! PSolver will overwrite rhopot_data with the potential.
    ! dummy_array for pot_ion (not used here as ixc=0 and sumpion=false)
    ! For actual use, pot_ion might need proper allocation if ixc/=0
    call PSolver(geometry_code, 'G', my_rank, num_procs, &
                 nx_dim, ny_dim, nz_dim, xc_functional_code, &
                 hx_spacing, hy_spacing, hz_spacing, &
                 rhopot_data, kernel_array, rhopot_data(1,1,1) , & ! Dummy pot_ion
                 hartree_energy, dummy_exc, dummy_vxc, &
                 offset_val, sum_ionic_potential, spin_components)

    if (my_rank == 0) then
      print *, "Poisson solve complete."
      print *, "Hartree Energy: ", hartree_energy
      ! rhopot_data now contains the electrostatic potential.
      ! ... (Process the potential data) ...
    end if

    ! Deallocate kernel
    if (associated(kernel_array)) deallocate(kernel_array)
    deallocate(rhopot_data)

    call MPI_Finalize(ierr)

  end subroutine solve_density_for_potential

end module MyPoissonProblem

program TestPoisson
  use MyPoissonProblem
  implicit none
  call solve_density_for_potential()
end program TestPoisson
```
