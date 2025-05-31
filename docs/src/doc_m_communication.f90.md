# `m_communication.f90`

## Overview

The `m_communication` module is central to enabling parallel execution of the `m_octree_mg` library. It encapsulates the necessary Message Passing Interface (MPI) functionalities to facilitate data exchange between different MPI processes. This module provides the mechanisms for common parallel tasks such as updating ghost cell regions (halo data) around local subdomains and transferring data segments during distributed multigrid operations like restriction and prolongation. While this module provides the generic transfer capability, specific data packing and unpacking for such operations are typically handled by other modules that utilize these communication routines.

## Key Components

### Modules

- **`m_communication`:** This module contains the subroutines and MPI interactions for inter-process data transfer.

### Functions/Subroutines

- **`mg_comm_init(mg, comm)`:**
  - **Description:** Initializes the MPI environment if it has not been initialized already. It then populates the `mg_t` data structure with essential MPI information: the MPI communicator (`mg%comm`), the rank of the current process (`mg%my_rank`), and the total number of processes (`mg%n_cpu`) in the communicator. An optional `comm` argument allows the user to specify a particular MPI communicator; if omitted, `MPI_COMM_WORLD` is used by default.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The main multigrid data structure where MPI information will be stored.
    - `comm (integer, intent(in), optional)`: An existing MPI communicator to be used. Defaults to `MPI_COMM_WORLD`.

- **`sort_and_transfer_buffers(mg, dsize)`:**
  - **Description:** This is the core routine for performing bulk data exchange between MPI processes. For each other process `i`:
    1. If data is queued for sending to process `i` (i.e., `mg%buf(i)%i_send > 0`), the send buffer `mg%buf(i)%send` is first sorted by the `sort_sendbuf` routine based on `mg%buf(i)%ix`. An asynchronous MPI send (`MPI_ISEND`) is then posted.
    2. If data is expected from process `i` (i.e., `mg%buf(i)%i_recv > 0`), an asynchronous MPI receive (`MPI_IRECV`) is posted for `mg%buf(i)%recv`.
    After initiating all non-blocking send and receive operations, the subroutine waits for all receive operations to complete, followed by waiting for all send operations to complete, using `MPI_WAITALL`. The MPI tag for these operations is hardcoded to `0`.
  - **Arguments:**
    - `mg (type(mg_t), intent(inout))`: The multigrid data structure containing the communication buffers (`mg%buf`).
    - `dsize (integer, intent(in))`: The size (number of `real(dp)` elements) of each logical item within the send/receive buffers. This is used by `sort_sendbuf` to correctly permute data.

### Private Helper Subroutines

- **`sort_sendbuf(gc, dsize)`:**
  - **Description:** Sorts the data in the send buffer `gc%send` according to the order specified by the `gc%ix` array (which contains original indices or keys for sorting). This is crucial if data items were added to the send buffer in an arbitrary order but need to be received in a specific, agreed-upon order. The sorting is performed out-of-place using a copy, and the `mrgrnk` routine (presumably a merge-rank sort algorithm from module `m_mrgrnk`) is used to obtain the permutation indices for reordering both `gc%send` and `gc%ix`.
  - **Arguments:**
    - `gc (type(mg_buf_t), intent(inout))`: The communication buffer structure for a specific partner rank, containing `send`, `ix`, and `i_send`, `i_ix` counters.
    - `dsize (integer, intent(in))`: The size of each data item in `gc%send` in terms of `real(dp)` elements.

## Important Variables/Constants

- **From `mg_t` (populated by `mg_comm_init` or used by routines):**
  - `mg%comm (integer)`: Stores the MPI communicator handle.
  - `mg%my_rank (integer)`: Rank of the current MPI process.
  - `mg%n_cpu (integer)`: Total number of MPI processes.
  - `mg%buf (array of type(mg_buf_t))`: Array of buffer structures, one for each MPI rank, used for sending and receiving data.
- **Standard MPI Constants (used from `mpi` module):**
  - `MPI_COMM_WORLD`: Default MPI communicator.
  - `MPI_DOUBLE`: MPI datatype for `real(dp)` variables.
  - `MPI_STATUSES_IGNORE`: Used in `MPI_WAITALL` when detailed status of each operation is not needed.
  - MPI error codes (`ierr`).
- **Communication Tag:** A default tag of `0` is used for messages in `sort_and_transfer_buffers`.

## Usage Examples

```fortran
! Conceptual example: Initialize MPI and perform a generic data exchange.

use m_data_structures
use m_communication
! Assuming 'mpi' module is available for MPI constants if not re-exported.

type(mg_t) :: my_multigrid_setup
integer    :: items_per_buffer_element

! 1. Initialize communication aspects of the multigrid setup
! This is typically done once at the beginning of a parallel run.
call mg_comm_init(my_multigrid_setup) ! Uses MPI_COMM_WORLD by default.
                                     ! Alternatively, pass a specific communicator:
                                     ! call mg_comm_init(my_multigrid_setup, my_custom_comm)

! Print rank for verification
print *, "Process ", my_multigrid_setup%my_rank, " of ", my_multigrid_setup%n_cpu, " initialized."

if (my_multigrid_setup%n_cpu > 1) then
    ! 2. Prepare data for sending and specify expected receive counts.
    ! This part is typically handled by other modules (e.g., m_ghost_cells, m_restrict).
    ! For example, a hypothetical routine in another module might:
    ! call prepare_data_for_exchange(my_multigrid_setup, items_per_buffer_element)
    ! This routine would:
    !   - Loop through destination ranks `dest_rank = 0, my_multigrid_setup%n_cpu - 1`.
    !   - Fill `my_multigrid_setup%buf(dest_rank)%send(1:count_send)` with data.
    !   - Fill `my_multigrid_setup%buf(dest_rank)%ix(1:count_ix)` with sorting keys.
    !   - Set `my_multigrid_setup%buf(dest_rank)%i_send = count_send`.
    !   - Set `my_multigrid_setup%buf(dest_rank)%i_ix = count_ix`.
    !   - Set `my_multigrid_setup%buf(source_rank)%i_recv = count_recv` for expected data.

    items_per_buffer_element = 1 ! Or mg_num_vars, etc., depending on what's being sent

    ! 3. Execute the sorted data transfer.
    ! Data in send buffers is sorted, then exchanged using MPI_Isend/Irecv, followed by MPI_Waitall.
    call sort_and_transfer_buffers(my_multigrid_setup, items_per_buffer_element)

    ! 4. Process the received data.
    ! Again, this would be handled by the modules that initiated the exchange.
    ! For example:
    ! call unpack_and_use_exchanged_data(my_multigrid_setup, items_per_buffer_element)
    ! This routine would read from `my_multigrid_setup%buf(source_rank)%recv(1:count_recv)`.

    print *, "Process ", my_multigrid_setup%my_rank, ": Buffer transfer complete."
else
    print *, "Single process run, no MPI communication performed by sort_and_transfer_buffers."
end if

! MPI_Finalize would be called elsewhere, typically at the very end of the program.
```

## Dependencies and Interactions

- **Internal Dependencies:**
  - **`m_data_structures`:** Absolutely essential. This module defines `type(mg_t)` (which stores MPI handles like `comm`, `my_rank`, `n_cpu`, and the array of communication buffers `buf`) and `type(mg_buf_t)` (which defines the structure of individual send/receive buffers: `send`, `recv`, `ix` arrays, and their respective counters `i_send`, `i_recv`, `i_ix`).
  - **`m_mrgrnk` (via `sort_sendbuf`):** This module provides the `mrgrnk` subroutine, which is used to perform a sort (likely a merge-rank sort) of indices to reorder the send buffer before transmission. This ensures data arrives in the order expected by the receiver.
- **External Libraries:**
  - **`mpi`:** The entire module is a wrapper and utility layer over the MPI library. All core communication operations (`MPI_Init`, `MPI_Comm_rank`, `MPI_Isend`, `MPI_Irecv`, `MPI_Waitall`, etc.) are direct calls to the MPI library.
- **Interactions with Other Components:**
  - **Program Initialization:** `mg_comm_init` is called during the initial setup phase of the simulation to prepare the MPI environment and store necessary MPI-related information in the `mg_t` structure.
  - **Data Preparation Modules (e.g., `m_ghost_cells`, `m_restrict`, `m_prolong`):** These modules are responsible for identifying what data needs to be exchanged. They would:
    1.  Populate the `send` array and `ix` array within `mg%buf(target_rank)` for each process that needs to receive data.
    2.  Set the `i_send` and `i_ix` counters for outgoing messages.
    3.  Set the `i_recv` counter for expected incoming messages.
    Once data is packed into buffers, these modules would then call `sort_and_transfer_buffers` to execute the physical data transfer.
  - **Data Consumption Modules:** After `sort_and_transfer_buffers` returns, the same modules (or related ones) would unpack the data from the `recv` array of `mg%buf(source_rank)` and use it for their intended purpose (e.g., updating ghost cell values in the local grid, providing data for restriction/prolongation calculations).
  - **Numerical Kernels (e.g., in `m_ahelmholtz`, `m_poisson`):** Any stencil-based computations or other operations that depend on data from neighboring domains rely on the ghost cells being correctly and timely updated. The communication facilitated by this module is thus a prerequisite for the correctness of these numerical operations in a parallel setting.
```
