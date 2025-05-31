# `psolver_main.f90`

## Overview

The `psolver_main.f90` file is a central piece of the `Poisson_Solver` module, designed to be `include`d within the main `poisson_solver.f90` file. Its primary contribution is the `PSolver` subroutine, which orchestrates the entire process of solving Poisson's equation ($\nabla^2 V = -4\pi\rho$, or variations depending on internal scaling) using Fast Fourier Transform (FFT) techniques.

`PSolver` is a versatile routine that handles various boundary conditions (Periodic, Surface, Free), different data distribution schemes for parallel processing (Global or Distributed), and can optionally incorporate exchange-correlation (XC) potentials and energies by interfacing with external XC libraries (like ABINIT's). This file also contains crucial helper subroutines for determining the appropriate dimensions for FFT grids and data arrays based on the problem's configuration: `PS_dim4allocation`, and the boundary-specific `P_FFT_dimensions`, `S_FFT_dimensions`, and `F_FFT_dimensions`.

## Key Components

### File Type
A collection of Fortran subroutines, intended for direct inclusion into `poisson_solver.f90`. It does not define a standalone Fortran `module`.

### Main Subroutine

-   **`PSolver(geocode, datacode, iproc, nproc, n01, n02, n03, ixc, hx, hy, hz, rhopot, karray, pot_ion, eh, exc, vxc, offset, sumpion, nspin)`**:
    -   **Description:** This is the main driver routine for solving Poisson's equation. Its operation can be summarized as:
        1.  **Dimension Calculation:** Determines optimal FFT grid dimensions (`n1,n2,n3`), padded real-space dimensions (`md1,md2,md3`), and kernel dimensions (`nd1,nd2,nd3`) by calling the appropriate `P/S/F_FFT_dimensions` routine based on the `geocode`.
        2.  **Memory Allocation:** Allocates local work arrays, primarily `zf` (to hold the density for FFT and then the potential) and `zfionxc` (for XC/ionic potentials).
        3.  **Exchange-Correlation (Optional):** If `ixc /= 0`, it calls `xc_energy` (from `xcenergy.f90`). `xc_energy` computes the XC energy and potential, potentially involving gradient calculations (via `3dgradient.f90`) if a GGA functional is selected. The input density from `rhopot` is appropriately processed (e.g., distributed, padded) and placed into `zf`.
        4.  **Core FFT Solve:** Calls the boundary-condition-specific solver routine (`P_PoissonSolver`, `S_PoissonSolver`, or `F_PoissonSolver` from `psolver_base.f90`). These routines perform the three main steps:
            a.  Forward 3D FFT of the density (from `zf`).
            b.  Multiplication in reciprocal space by the pre-computed `karray` (the $G(k)$ kernel).
            c.  Inverse 3D FFT to transform the result back to real space, yielding the Hartree potential (stored back in `zf`).
        5.  **Potential Correction & Combination:** The computed Hartree potential is corrected (e.g., for `offset` in periodic systems). If XC calculations were performed, the XC potential (from `zfionxc`) is added. If `sumpion` is true, the external `pot_ion` is also added. The final combined potential is stored in the output `rhopot` array.
        6.  **Energy Calculation:** The Hartree energy (`eh`) and XC related energies (`exc`, `vxc`) are computed locally and then summed across all MPI processes using `MPI_ALLREDUCE`.
        7.  **Data Gathering (Optional):** If `datacode = 'G'` (global data), the final potential in `rhopot` (and `pot_ion` if modified) is gathered from all MPI processes to all processes using `MPI_ALLGATHERV`.
    -   **Key Arguments:**
        -   `geocode (character(len=1), intent(in))`: Boundary condition type ('F', 'S', 'P').
        -   `datacode (character(len=1), intent(in))`: Data distribution type ('G' for global, 'D' for distributed).
        -   `iproc, nproc (integer, intent(in))`: MPI rank of the current process and total number of processes.
        -   `n01, n02, n03 (integer, intent(in))`: Global dimensions of the original real-space grid.
        -   `ixc (integer, intent(in))`: Code for the exchange-correlation functional (0 for none; ABINIT convention).
        -   `hx, hy, hz (real(kind=8), intent(in))`: Grid spacings in x, y, z directions.
        -   `rhopot (real(kind=8), dimension(*), intent(inout))`: On input, this array holds the source density $\rho$. On output, it is overwritten with the computed potential $V_H$ (plus $V_{XC}$ and $V_{ion}$ if applicable). Its actual dimensions depend on `datacode` and must be consistent with `PS_dim4allocation`.
        -   `karray (real(kind=8), dimension(*), intent(in))`: The pre-computed Poisson kernel $G(k)$ in reciprocal space, obtained from `createKernel`.
        -   `pot_ion (real(kind=8), dimension(*), intent(inout))`: If `ixc /= 0` and `sumpion` is `.false.`, this array is overwritten with $V_{XC}$. If `sumpion` is `.true.`, this array provides an input ionic/external potential that is added to the final result in `rhopot`.
        -   `eh (real(kind=8), intent(out))`: Calculated Hartree energy ($1/2 \int \rho V_H d\mathbf{x}$).
        -   `exc (real(kind=8), intent(out))`: Calculated exchange-correlation energy $E_{XC}$.
        -   `vxc (real(kind=8), intent(out))`: Calculated $\int \rho V_{XC} d\mathbf{x}$.
        -   `offset (real(kind=8), intent(in))`: A constant potential offset to be applied, primarily for periodic systems.
        -   `sumpion (logical, intent(in))`: If `.true.` and `ixc /=0`, `pot_ion` (input) is added to the final potential in `rhopot`.
        -   `nspin (integer, intent(in))`: Number of spin components (1 for non-spin-polarized, 2 for spin-polarized).

### Dimension Calculation Subroutines

-   **`PS_dim4allocation(geocode, datacode, iproc, nproc, n01, n02, n03, ixc, n3d, n3p, n3pi, i3xcsh, i3s)`**:
    -   **Description:** This routine is crucial for users to determine the correct sizes for allocating the `rhopot` and `pot_ion` arrays, especially for distributed data (`datacode='D'`) and when GGA XC functionals (`ixc >= 11`) are used, as these require larger halo regions for gradient computations. It calls one of `P/S/F_FFT_dimensions` based on `geocode` to get base FFT and real-space dimensions, then adjusts these based on `ixc` and `datacode`.
    -   **Key Output Arguments:**
        -   `n3d (integer, intent(out))`: The required size of the third dimension for the `rhopot` array (holding density) on the current process.
        -   `n3p (integer, intent(out))`: The required size of the third dimension for the `rhopot` array (when it holds potential) on the current process.
        -   `n3pi (integer, intent(out))`: The required size of the third dimension for the `pot_ion` array on the current process.
        -   `i3xcsh (integer, intent(out))`: The shift in the third dimension index needed to access the non-overlapping part of the potential when XC gradients are involved.
        -   `i3s (integer, intent(out))`: The starting global plane index in the third dimension for the data slice handled by the current process.

-   **`P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)`**
-   **`S_FFT_dimensions(...)`**
-   **`F_FFT_dimensions(...)`**
    -   **Description:** These helper subroutines calculate various sets of dimensions required for FFTs under Periodic, Surface, or Free boundary conditions, respectively.
        -   `m1,m2,m3`: Real-space dimensions, with `m2=n03` and `m3=n02` (a convention used internally, possibly y-z swap).
        -   `n1,n2,n3`: Padded dimensions suitable for FFT (e.g., products of small primes from `fourier_dim`). For 'F' and 'S', these are typically doubled for zero-padding or handling non-periodic directions.
        -   `md1,md2,md3`: Dimensions of the real-space data arrays used in FFTs, potentially adjusted for MPI distribution (e.g., `md2` might be made a multiple of `nproc`).
        -   `nd1,nd2,nd3`: Dimensions for the kernel array in reciprocal space (often half-dimensions due to symmetries), adjusted for MPI distribution.
    -   They utilize `fourier_dim` (from `fft3d.f90`) to find FFT-compatible lengths.

## Important Variables/Constants

-   **`nordgr (integer, parameter)` in `PSolver`**: Defines the order for finite-difference gradient calculations (fixed to 4), relevant for GGA XC functionals.
-   **Internal work arrays `zf` and `zfionxc`**: Dynamically allocated within `PSolver` to hold data during FFT operations and for XC potential components.
-   **MPI tags**: Standard MPI tag 0 is used for communications within `PSolver`.

## Dependencies and Interactions

-   **Internal (within `Poisson_Solver` module via `include` mechanism):**
    -   `PSolver` calls `P/S/F_FFT_dimensions` (defined in the same file).
    -   If `ixc /= 0`, `PSolver` calls `xc_energy` (from `xcenergy.f90`).
    -   `PSolver` calls one of `P_PoissonSolver`, `S_PoissonSolver`, or `F_PoissonSolver` (from `psolver_base.f90`) to perform the core FFT-based solution.
    -   The `P/S/F_FFT_dimensions` routines call `fourier_dim` (from `fft3d.f90`).
-   **External Libraries:**
    -   **MPI:** Heavily used in `PSolver` for global communication: `MPI_BCAST` (for potential offset in periodic case), `MPI_ALLREDUCE` (for summing energies), and `MPI_ALLGATHERV` (if `datacode='G'` to ensure all processes have the full result).
-   **Interactions with Other Components:**
    -   `PSolver` is the primary computational routine of the `Poisson_Solver` module. It expects the $G(k)$ kernel to have been pre-calculated by `createKernel` (from `build_kernel.f90`).
    -   The dimensions of input arrays `rhopot` and `pot_ion` must be determined using `PS_dim4allocation` to ensure compatibility, especially in parallel and with XC calculations.
    -   This entire `Poisson_Solver` package is used by the `m_free_space` module within the `m_octree_mg` library to compute the far-field potential for open boundary conditions, typically on the coarsest grid of the multigrid hierarchy.

## Usage Examples

The `PSolver` subroutine is complex due to its many options. A typical invocation sequence involves first determining array sizes, then creating the kernel, then calling the solver.

```fortran
! Conceptual Example: Using PSolver after kernel creation.

module ExampleFullPoissonSolve
  use Poisson_Solver ! Provides PSolver, PS_dim4allocation, createKernel
  use mpi
  implicit none

  subroutine run_fft_poisson_solve()
    real(kind=8), pointer :: kernel_data(:)
    real(kind=8), dimension(:,:,:), allocatable :: density_potential_3d_array
    real(kind=8), dimension(:,:,:), allocatable :: ionic_potential_3d_array ! Or for Vxc output
    real(kind=8) :: hx, hy, hz, hartree_energy, xc_energy_val, vxc_integral_val, potential_offset
    integer :: nx, ny, nz, mpi_rank_id, num_mpi_procs, xc_func_id
    integer :: n3d_local, n3p_local, n3pi_local, i3xc_shift, i3_start_plane, mpi_err
    character(len=1) :: bc_geocode, data_layout_code
    logical :: add_ionic_potential_to_hartree

    ! 1. Initialize MPI and basic parameters
    call MPI_Init(mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank_id, mpi_err)
    call MPI_Comm_size(MPI_COMM_WORLD, num_mpi_procs, mpi_err)

    nx = 96; ny = 96; nz = 96       ! Global grid dimensions
    hx = 0.2d0; hy = 0.2d0; hz = 0.2d0 ! Grid spacings
    bc_geocode = 'F'              ! Free boundary conditions
    data_layout_code = 'D'        ! Distributed data
    xc_func_id = 0                ! 0 for no exchange-correlation
    add_ionic_potential_to_hartree = .false.
    potential_offset = 0.0d0      ! Usually for periodic systems

    ! 2. Determine allocation sizes for distributed arrays
    call PS_dim4allocation(bc_geocode, data_layout_code, mpi_rank_id, num_mpi_procs, &
                           nx, ny, nz, xc_func_id, &
                           n3d_local, n3p_local, n3pi_local, i3xc_shift, i3_start_plane)

    ! 3. Allocate arrays based on determined local sizes
    !    Note: Dimensions might be swapped internally by FFT routines (n02, n03 often swapped)
    !    For simplicity, using n01, n02, n3d_local (actual mapping is complex)
    allocate(density_potential_3d_array(nx, ny, n3d_local)) ! Using n01,n02 for first two dims
    allocate(ionic_potential_3d_array(nx, ny, n3pi_local))  ! Same for pot_ion
    density_potential_3d_array = 0.0d0
    ionic_potential_3d_array = 0.0d0
    ! ... Fill the local slice of density_potential_3d_array with density values ...
    ! ... (from i3_start_plane to i3_start_plane + n3d_local - 1) ...

    ! 4. Create the Poisson kernel
    !    (itype_scf=14 is a common choice for scaling functions)
    call createKernel(bc_geocode, nx, ny, nz, hx, hy, hz, &
                      14, mpi_rank_id, num_mpi_procs, kernel_data)

    ! 5. Call the Poisson Solver
    call PSolver(bc_geocode, data_layout_code, mpi_rank_id, num_mpi_procs, &
                 nx, ny, nz, xc_func_id, hx, hy, hz, &
                 density_potential_3d_array, kernel_data, ionic_potential_3d_array, &
                 hartree_energy, xc_energy_val, vxc_integral_val, &
                 potential_offset, add_ionic_potential_to_hartree, 1) ! nspin=1

    ! 6. Results
    !    density_potential_3d_array now contains the computed electrostatic potential
    !    (potentially plus XC and ionic contributions based on flags).
    !    hartree_energy, xc_energy_val, vxc_integral_val contain calculated energies.
    if (mpi_rank_id == 0) then
      print *, "PSolver finished. Hartree Energy:", hartree_energy
    end if

    ! 7. Deallocate
    if (associated(kernel_data)) deallocate(kernel_data)
    deallocate(density_potential_3d_array)
    deallocate(ionic_potential_3d_array)

    call MPI_Finalize(mpi_err)

  end subroutine run_fft_poisson_solve
end module ExampleFullPoissonSolve

program TestPoissonSolverMain
  use ExampleFullPoissonSolve
  implicit none
  call run_fft_poisson_solve()
end program TestPoissonSolverMain
```
