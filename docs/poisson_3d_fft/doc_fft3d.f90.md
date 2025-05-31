# `fft3d.f90`

## Overview

The file `fft3d.f90` is a collection of low-level Fortran subroutines that provide the core Fast Fourier Transform (FFT) capabilities used within the `Poisson_Solver` module of the `poisson_3d_fft` package. This file does not define a self-contained Fortran `module` in the typical sense; rather, its subroutines are designed to be directly `include`d into other source files where FFT functionality is required (e.g., `psolver_base.f90`, `build_kernel.f90`).

The routines handle critical aspects of the FFT process:
1.  Determining optimal grid dimensions that are efficient for FFT algorithms (products of small prime factors).
2.  Pre-calculating the trigonometric coefficients (twiddle factors) needed for the FFT stages.
3.  Executing the individual 1D complex-to-complex FFT steps (passes) that form the basis of multi-dimensional FFTs.

This appears to be a custom implementation of the Cooley-Tukey FFT algorithm, or a similar radix-based approach, rather than a wrapper around a standard library like FFTW.

## Key Components

### File Type
A collection of Fortran subroutines, intended for direct inclusion. No `module` statement is present in this file.

### Key Subroutines (Effectively Public via `include`)

-   **`fourier_dim(n, n_next)`**:
    -   **Description:** Given an input integer `n`, this subroutine finds and returns in `n_next` the smallest integer greater than or equal to `n` that is suitable for efficient FFT computation. It uses a hardcoded list (`idata`) of integers that are products of small prime factors (2, 3, 5), which are known to be efficient lengths for FFT algorithms.
    -   **Arguments:**
        -   `n (integer, intent(in))`: The desired minimum dimension.
        -   `n_next (integer, intent(out))`: The returned FFT-friendly dimension.

-   **`ctrig(n, trig, after, before, now, isign, ic)`**:
    -   **Description:** This subroutine prepares the necessary trigonometric coefficients (twiddle factors, which are complex exponentials $e^{\pm 2\pi i k/N}$) for performing an FFT of length `n`. It also determines the factorization of `n` into stages suitable for the Cooley-Tukey algorithm. The factorization information (radices for each stage, number of transforms at each stage, etc.) is stored in the `after`, `before`, and `now` arrays, with `ic` being the number of stages. The `isign` argument (+1 or -1) determines the sign in the exponent, allowing for both forward and inverse transforms. The pre-calculated sine and cosine values are stored in the `trig` array.
    -   **Arguments:**
        -   `n (integer)`: The length of the 1D transform.
        -   `trig (real(kind=8), dimension(2,nfft_max))`: Output array storing the real and imaginary parts of the twiddle factors.
        -   `after(7), before(7), now(7) (integer)`: Output arrays describing the factorization of `n` for the FFT stages.
        -   `isign (integer)`: Sign for the exponent in twiddle factors (e.g., -1 for forward, +1 for inverse, or vice-versa depending on convention).
        -   `ic (integer)`: Output, the number of FFT stages.

-   **`fftstp(mm, nfft, m, nn, n, zin, zout, trig, after, now, before, isign)`**:
    -   **Description:** This subroutine performs a single "step" or "pass" of a 1D complex-to-complex FFT. It processes `nfft` transforms simultaneously. The `now` argument (from `ctrig`) specifies the radix of the current FFT stage (e.g., 2, 3, 4, 5, 6, 8). The routine uses the pre-calculated twiddle factors from `trig` and the stage information (`after`, `before`) to combine input data from `zin` and store the result in `zout`. This routine is called multiple times (once for each stage `ic` determined by `ctrig`) to complete a full 1D FFT. The `isign` parameter ensures the correct twiddle factors are used for either a forward or an inverse transform.
    -   **Arguments:**
        -   `mm, nfft, m, nn, n (integer)`: Parameters related to data layout, number of transforms, and dimensions for efficient processing (often related to cache blocking).
        -   `zin (real(kind=8), dimension(2,mm,m))`: Input complex data (real and imaginary parts).
        -   `zout (real(kind=8), dimension(2,nn,n))`: Output complex data.
        -   `trig (real(kind=8), dimension(2,nfft_max))`: Pre-calculated twiddle factors from `ctrig`.
        -   `after, now, before (integer)`: Current FFT stage parameters from `ctrig`.
        -   `isign (integer)`: Sign for the transform direction.

## Important Variables/Constants

-   **`nfft_max (integer, parameter)`**: Defined in `ctrig` and `fftstp` (value 24000). This limits the maximum length of a 1D FFT for which twiddle factors can be pre-stored in the `trig` array.
-   **`idata (integer, dimension(ndata), parameter)`**: In `fourier_dim`, a hardcoded list of 180 integers that are products of small primes (2, 3, 5), suitable for FFTs. The largest is 24000.
-   **`idata (integer, dimension(7,ndata), parameter)`**: In `ctrig`, a hardcoded table providing factorizations for the integers listed in `fourier_dim`'s `idata`. These factorizations (into factors of 2, 3, 4, 5, 6, 8) determine the stages of the `fftstp` routine.
-   **Numerical Constants:** `rt2i = sqrt(0.5d0)` (i.e., $1/\sqrt{2}$), `cos2/sin2` (for $\cos(2\pi/5), \sin(2\pi/5)$), `cos4/sin4` (for $\cos(4\pi/5), \sin(4\pi/5)$), `bb` (for $\sin(2\pi/3)$ or $\cos(\pi/6)$ related to factor-3 FFTs) are used in `fftstp` for specific radices.

## Dependencies and Interactions

-   **Internal Dependencies:**
    -   `fftstp` relies on the setup performed by `ctrig` (factorization of N and twiddle factors).
-   **External Libraries:**
    -   None. The FFT implementation appears to be entirely self-contained within this file and the routines that call it. It does not seem to wrap standard FFT libraries like FFTW.
-   **Interactions with Other Components (within `poisson_3d_fft`):**
    -   **`psolver_base.f90` (e.g., `P_PoissonSolver`, `S_PoissonSolver`, `F_PoissonSolver`) and `build_kernel.f90` (e.g., `kernelfft`):** These higher-level routines are the primary callers of `ctrig` and `fftstp`. They orchestrate the full 3D FFT by:
        1.  Calling `ctrig` to prepare for 1D transforms along each dimension.
        2.  Calling `fftstp` multiple times for each dimension to perform the 1D FFTs. This often involves data reordering or transposition steps (handled within `psolver_base.f90` or `kernelfft.f90` through helper routines like `switch`, `mpiswitch`, `scramble_unpack`, etc.) to present data to `fftstp` as contiguous 1D arrays or slices.
    -   **`psolver_main.f90` (via `P/S/F_FFT_dimensions`):** Uses `fourier_dim` to ensure that the grid dimensions chosen for the Poisson solve are efficient for the FFT routines provided in this file.

## Usage Examples

The subroutines in `fft3d.f90` are low-level components of the 3D FFT algorithm within the `Poisson_Solver` module. A typical end-user of the `Poisson_Solver` (or the `m_octree_mg` library) would not call these routines directly. They are used internally by routines like `P_PoissonSolver`, `S_PoissonSolver`, `F_PoissonSolver`, and `kernelfft`.

**Conceptual Invocation for a 1D FFT (as used by higher-level 3D FFT routines):**

```fortran
! Conceptual example of how a higher-level routine might use ctrig and fftstp
! to perform a 1D FFT on multiple data vectors.

implicit none
real(kind=8) :: my_data_array(2, 128, 10) ! Example: 10 complex vectors, each of length 128
real(kind=8) :: work_array(2, 128, 10)    ! Work array for ping-ponging
real(kind=8) :: twiddle_factors(2, 24000) ! Max size from nfft_max
integer :: factors_after(7), factors_now(7), factors_before(7)
integer :: num_fft_stages, current_stage
integer :: N_length, num_vectors
integer :: isign_direction

N_length = 128
num_vectors = 10 ! This would be 'nfft' in fftstp if processing all vectors in one call
isign_direction = -1 ! For a forward FFT (convention-dependent)

! 1. Pre-calculate twiddle factors and FFT stages for length N_length
call ctrig(N_length, twiddle_factors, factors_after, factors_before, &
           factors_now, isign_direction, num_fft_stages)

! 2. Perform the FFT using fftstp for each stage
!    Source and destination arrays might ping-pong between my_data_array and work_array.
!    For simplicity, assume source is my_data_array, destination is work_array for 1st stage.

! Loop over each stage of the FFT
do current_stage = 1, num_fft_stages
  if (mod(current_stage, 2) == 1) then ! Odd stages: my_data_array -> work_array
    call fftstp(N_length, num_vectors, factors_now(current_stage), & ! mm, nfft, m (radix)
                N_length, num_vectors,                               & ! nn, n
                my_data_array, work_array,                            & ! zin, zout
                twiddle_factors,                                      & ! trig
                factors_after(current_stage), factors_now(current_stage), &
                factors_before(current_stage), isign_direction)
  else ! Even stages: work_array -> my_data_array
    call fftstp(N_length, num_vectors, factors_now(current_stage), &
                N_length, num_vectors,                               &
                work_array, my_data_array,                            &
                twiddle_factors,                                      &
                factors_after(current_stage), factors_now(current_stage), &
                factors_before(current_stage), isign_direction)
  end if
end do

! Final result is in my_data_array if num_fft_stages is even,
! or in work_array if num_fft_stages is odd. A final copy might be needed.

! To find an FFT-friendly dimension:
integer :: my_dim, good_fft_dim
my_dim = 120
call fourier_dim(my_dim, good_fft_dim)
! good_fft_dim will be 120

my_dim = 121
call fourier_dim(my_dim, good_fft_dim)
! good_fft_dim will be 125
```
