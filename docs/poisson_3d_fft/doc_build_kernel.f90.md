# `build_kernel.f90`

## Overview

The `build_kernel.f90` file is a crucial part of the `Poisson_Solver` module within the `poisson_3d_fft` package. It is not a standalone Fortran module but a collection of subroutines intended to be directly `include`d into the main `Poisson_Solver` module file (`poisson_solver.f90`). Its primary responsibility is to construct the Green's function (often referred to as the Poisson kernel) in reciprocal (Fourier) space, denoted as $G(k)$. This kernel is essential for solving Poisson's equation $\nabla^2 V = -4\pi\rho$ via FFT methods, where the solution in Fourier space is $V(k) = G(k) \rho(k)$.

The specific form of $G(k)$ and the method of its construction depend significantly on the boundary conditions of the problem ('F' for Free, 'S' for Surface, 'P' for Periodic) and, for non-periodic cases, on the type of interpolating scaling functions used to represent the Green's function.

## Key Components

### File Type
A collection of Fortran subroutines, designed for direct inclusion into `poisson_solver.f90`.

### Main Public-Facing Subroutine (Exposed via `Poisson_Solver`)

-   **`createKernel(geocode, n01, n02, n03, hx, hy, hz, itype_scf, iproc, nproc, kernel)`**:
    -   **Description:** This is the primary routine made available through the `Poisson_Solver` module that originates from this file. It acts as a dispatcher based on the `geocode` argument:
        -   If `geocode = 'P'` (Periodic): The kernel is simple ($1/k^2$, often handled directly in the solver), so this routine may only allocate a minimal placeholder for `kernel` or perform basic setup. The actual $1/k^2$ scaling is applied in `P_PoissonSolver`.
        -   If `geocode = 'S'` (Surface): It calls `S_FFT_dimensions` (from `psolver_main.f90`) to determine appropriate FFT grid dimensions and then invokes the `Surfaces_Kernel` subroutine to compute the specific kernel $G(k)$ for surface boundary conditions.
        -   If `geocode = 'F'` (Free): It calls `F_FFT_dimensions` (from `psolver_main.f90`) and then invokes the `Free_Kernel` subroutine to compute $G(k)$ for free-space (isolated system) boundary conditions.
    -   **Output:** The `kernel` argument, a pointer, is allocated and filled with the computed reciprocal space Green's function values. For parallel runs (`nproc > 1`), this kernel data is typically distributed across processes.
    -   **Arguments:**
        -   `geocode (character(len=1), intent(in))`: Specifies boundary conditions ('P', 'S', 'F').
        -   `n01, n02, n03 (integer, intent(in))`: Dimensions of the real-space grid.
        -   `hx, hy, hz (real(kind=8), intent(in))`: Grid spacings.
        -   `itype_scf (integer, intent(in))`: Order of interpolating scaling functions used for kernel construction in non-periodic cases.
        -   `iproc, nproc (integer, intent(in))`: MPI rank of the current process and total number of processes.
        -   `kernel (real(kind=8), pointer, intent(out))`: On exit, points to the allocated array containing the computed kernel $G(k)$.

### Core Kernel Construction Subroutines (Internal to `build_kernel.f90`)

-   **`Surfaces_Kernel(n1, n2, n3, m3, nker1, nker2, nker3, h1, h2, h3, itype_scf, karray, iproc, nproc)`**:
    -   **Description:** Computes the Poisson kernel for surface boundary conditions (periodic in xz, isolated in y, where y is the third dimension `n3`/`m3` due to internal dimension swapping). This involves:
        1.  Generating 1D Green's function components in the non-periodic direction using interpolating scaling functions (via `calculates_green_opt_muzero` for $k_{y_p}=0$ and `calculates_green_opt` for $k_{y_p} \neq 0$, where $k_{y_p}$ is the wavevector component in the periodic plane). These helpers utilize `scaling_function` (from `scaling_function.f90`) and numerical integration with pre-defined polynomial coefficients (`cpol`).
        2.  Performing a "half-FFT" like procedure along the non-periodic dimension for each $(k_x, k_z)$ pair.
        3.  Distributing the final `karray` among MPI processes using `MPI_ALLTOALL`.
    -   **Key Techniques:** Use of scaling functions, specialized Green's function calculations for $k_p=0$ and $k_p \neq 0$, custom FFT-like processing for the mixed boundary condition.

-   **`Free_Kernel(n01, n02, n03, nfft1, nfft2, nfft3, n1k, n2k, n3k, hx, hy, hz, itype_scf, iproc, nproc, karray)`**:
    -   **Description:** Computes the Poisson kernel for free boundary conditions (isolated system). This is typically more complex than the periodic case. The method involves:
        1.  Representing the $1/r$ interaction using a sum of Gaussian functions via Gaussian quadrature (`gequad` provides abscissas `p_gauss` and weights `w_gauss`).
        2.  For each Gaussian, constructing its representation using interpolating scaling functions (via `scaling_function` and `scf_recursion` from `scaling_function.f90`). This effectively creates a real-space representation of the kernel (`kp`).
        3.  Transforming this real-space kernel `kp` to reciprocal space `karray` using a 3D FFT via the `kernelfft` subroutine.
    -   **Key Techniques:** Gaussian quadrature for $1/r$, representation with scaling functions, 3D FFT of the constructed real-space kernel.

### Helper Subroutines (Selected)

-   **`calculates_green_opt(...)`, `calculates_green_opt_muzero(...)`**: Used by `Surfaces_Kernel` to compute 1D Green's functions for the non-periodic direction using scaling functions and handling the $k_p=0$ (muzero) case separately.
-   **`gequad(nterms, p, w, urange, drange, acc)`**: Provides `nterms` abscissas (`p`) and weights (`w`) for Gaussian quadrature. Data for 89 points is hardcoded. Used in `Free_Kernel`.
-   **`kernelfft(n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,nproc,iproc,zf,zr)`**: Performs a 3D FFT. Takes a real-space kernel `zf` (e.g., `kp` from `Free_Kernel`) and computes its complex reciprocal-space counterpart `zr` (which becomes `karray`). It manages data distribution for MPI and calls 1D FFT routines (`fftstp` from `fft3d.f90`).
-   **`indices(...)`**: A small utility for calculating array indices, used within `Surfaces_Kernel`.
-   **`inserthalf(...)`, `realcopy(...)`, `switch(...)`, `mpiswitch(...)`**: Data manipulation routines used by `kernelfft` to arrange data for 1D FFT passes or MPI communication.

## Important Variables/Constants

-   **`n_points (integer, parameter)`**: Within `Surfaces_Kernel` and `Free_Kernel`, number of points for fine-grained integration related to scaling functions (value: $2^6=64$).
-   **`n_gauss (integer, parameter)`**: In `Free_Kernel`, the number of Gaussian functions used to represent $1/r$ (value: 89).
-   **`p_gauss, w_gauss (real(kind=8), dimension(n_gauss))`**: In `Free_Kernel`, arrays storing the exponents and weights for the Gaussian expansion of $1/r$, initialized by `gequad`.
-   **`cpol (real(kind=8), dimension(9,8))`**: In `Surfaces_Kernel`, hardcoded coefficients for polynomial interpolation used in the Green's function calculation helpers.
-   **`iproc_verbose (integer)`**: Module variable from `Poisson_Solver` used here to control print statements.

## Dependencies and Interactions

-   **Internal Dependencies (within `poisson_3d_fft` directory):**
    -   Relies on dimension calculation routines like `F_FFT_dimensions` and `S_FFT_dimensions` (defined in `psolver_main.f90` but available via include).
    -   `Surfaces_Kernel` and `Free_Kernel` depend on `scaling_function` and associated routines (like `scf_recursion_XX`) from `scaling_function.f90` for constructing parts of the Green's function.
    -   `Free_Kernel` (via `kernelfft`) and `Surfaces_Kernel` (for its "half-FFT") depend on the 1D FFT routines `fftstp` and `ctrig` from `fft3d.f90`.
-   **External Libraries:**
    -   **MPI:** The `MPI_ALLTOALL` routine is used in `Surfaces_Kernel` and `kernelfft` for distributing the computed kernel data among processes. `MPI_COMM_WORLD` is used.
-   **Interactions with Other Components:**
    -   The main entry point `createKernel` is called by the `Poisson_Solver` module (specifically, it's a public routine of `Poisson_Solver` whose implementation is in this file).
    -   The output `kernel` (pointer to `karray`) from `createKernel` is a critical input to the main `PSolver` subroutine (defined in `psolver_main.f90` and implemented in `psolver_base.f90`). `PSolver` uses this kernel to multiply the FFT of the density in reciprocal space.

## Usage Examples

The `createKernel` subroutine is typically not called directly by an end-user of the `m_octree_mg` library, but rather internally by `m_free_space` when it initializes the `Poisson_Solver`. However, a conceptual call would look like this:

```fortran
! Conceptual Example: How createKernel is invoked.

module ExampleKernelCreation
  use Poisson_Solver ! This makes createKernel (and its dependencies) available
  use mpi            ! MPI environment is needed by createKernel and its helpers
  implicit none

  subroutine generate_poisson_kernel()
    real(kind=8), pointer :: my_kernel_array(:)
    real(kind=8) :: hx, hy, hz       ! Grid spacings
    integer    :: nx, ny, nz       ! Grid dimensions
    integer    :: mpi_current_rank, num_mpi_tasks
    integer    :: scaling_func_order, mpi_error_code
    character(len=1) :: boundary_type_code

    ! Initialize MPI (usually done at a higher level)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_current_rank, mpi_error_code)
    call MPI_Comm_size(MPI_COMM_WORLD, num_mpi_tasks, mpi_error_code)

    ! Define problem parameters
    nx = 128; ny = 128; nz = 128
    hx = 0.1d0; hy = 0.1d0; hz = 0.1d0
    boundary_type_code = 'F'  ! Free-space boundary conditions
    scaling_func_order = 16   ! Order of interpolating scaling functions

    ! Call createKernel to compute and allocate the kernel
    call createKernel(geocode=boundary_type_code, &
                      n01=nx, n02=ny, n03=nz, &
                      hx=hx, hy=hy, hz=hz, &
                      itype_scf=scaling_func_order, &
                      iproc=mpi_current_rank, nproc=num_mpi_tasks, &
                      kernel=my_kernel_array)

    ! my_kernel_array is now allocated and contains the computed kernel in reciprocal space,
    ! distributed among MPI processes if num_mpi_tasks > 1.
    ! This kernel can then be passed to the PSolver routine.

    ! ... (Further operations, e.g., calling PSolver) ...

    ! Deallocate the kernel when it's no longer needed
    if (associated(my_kernel_array)) then
      deallocate(my_kernel_array)
    end if

  end subroutine generate_poisson_kernel
end module ExampleKernelCreation

program TestDriver
  use ExampleKernelCreation
  use mpi
  implicit none
  integer :: ierr
  call MPI_Init(ierr)
  call generate_poisson_kernel()
  call MPI_Finalize(ierr)
end program TestDriver
```
