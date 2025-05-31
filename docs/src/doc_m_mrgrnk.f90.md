# `m_mrgrnk.f90`

## Overview

The `m_mrgrnk` module provides a sorting utility. Contrary to what its name might initially suggest (e.g., related to MPI ranks), this module implements a "Merge-sort ranking" algorithm for an array of integers. Instead of returning the sorted data directly, it outputs an array of indices (ranks). These indices specify the order of the elements from the original input array if they were to be arranged in ascending sorted order. For example, the first element of the rank array gives the index of the smallest element in the input array, the second element of the rank array gives the index of the second smallest, and so on. This type of ranking is useful when the original positions of elements are important or when multiple arrays need to be sorted based on the keys in one primary array.

## Key Components

### Modules

- **`m_mrgrnk`:** Encapsulates the merge-sort ranking functionality.

### Public Interface/Subroutine

- **`mrgrnk(XDONT, IRNGT)`:**
  - **Description:** This is the sole public interface provided by the module, which internally calls the `I_mrgrnk` subroutine. It takes an array of integers and returns an array of indices that represent the sort order.
  - **Arguments:**
    - `XDONT (integer, dimension(:), intent(in))`: The input array of integer values to be ranked. The name `XDONT` (presumably "X data original not touched") is slightly unconventional for an input array that *is* the data to be sorted by rank.
    - `IRNGT (integer, dimension(:), intent(out))`: The output array, which will be populated with the 1-based indices of the `XDONT` array in ascending order of their values. For example, `IRNGT(1)` will hold the index `j` such that `XDONT(j)` is the smallest value in `XDONT`. The size of `IRNGT` processed will be `min(size(XDONT), size(IRNGT))`.

### Private Subroutine

- **`I_mrgrnk(XDONT, IRNGT)`:**
  - **Description:** This is the actual implementation of the merge-sort ranking algorithm. The algorithm works by iteratively merging ordered sub-sequences. The implementation includes specific optimizations for the initial passes (creating ordered couples, then ordered quadruplets) before entering a general loop that doubles the length of ordered subsets in each iteration. It uses an internal work array `JWRKT` (a copy of parts of `IRNGT` during merge steps) to facilitate the merge operations.
  - **Arguments:** Same as the public `mrgrnk` interface.

## Important Variables/Constants

- **`kdp (integer, parameter, private)`:** Defined as `selected_real_kind(15)`. This constant is declared in the module but is not used within the `I_mrgrnk` subroutine, as the sorting logic is purely for integer arrays. It might be a remnant from a previous version or intended for a different, unincluded routine.
- **`JWRKT (integer, dimension(SIZE(IRNGT)))` (local in `I_mrgrnk`):** A work array used during the merge phases of the sort to temporarily hold parts of the index array being merged.

## Dependencies and Interactions

- **Internal Dependencies:**
  - The module relies only on standard Fortran intrinsic functions like `min`, `size`, and `modulo`.
- **External Libraries:**
  - None. This is a self-contained sorting utility.
- **Interactions with Other Components:**
  - **`m_communication`:** The `sort_sendbuf` subroutine within `m_communication` uses `mrgrnk` to sort the `gc%ix` array (which likely contains indices or keys for data elements in a send buffer) and then reorders the actual send buffer `gc%send` according to these sorted ranks. This is essential for ensuring that data elements within MPI messages are in a consistent and expected order when received by other processes, especially if they were packed into the buffer in a different order.
  - **General Utility:** Any other module within the `m_octree_mg` library (or an external application using it) that needs to determine the sorted order of an integer array without modifying the array itself can use `m_mrgrnk`.

## Usage Examples

```fortran
! Example demonstrating the use of the mrgrnk subroutine.

use m_mrgrnk
implicit none

integer, dimension(7) :: my_data_array
integer, dimension(7) :: sorted_indices
integer :: i

! Initialize data
my_data_array = [45, 12, 89, 2, 67, 22, 50]

! Call mrgrnk to get the indices that would sort my_data_array
call mrgrnk(my_data_array, sorted_indices)

! Output the results
print "(A, 7I4)", "Original Data : ", my_data_array
print "(A, 7I4)", "Rank Indices  : ", sorted_indices

print *, ""
print *, "Data elements in sorted order (accessed via sorted_indices):"
do i = 1, size(my_data_array)
  print "(I4)", my_data_array(sorted_indices(i))
end do

! Expected output:
! Original Data :   45  12  89   2  67  22  50
! Rank Indices  :    4   2   6   1   7   5   3
!
! Data elements in sorted order (accessed via sorted_indices):
!   2
!  12
!  22
!  45
!  50
!  67
!  89
```
