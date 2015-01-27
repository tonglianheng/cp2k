MODULE qs_fb_atomic_matrix_methods

  USE kinds,                     ONLY: dp, &
                                       int_8
  USE qs_fb_env_types,           ONLY: fb_env_obj, &
                                       fb_env_get
  USE qs_fb_atomic_halo_types,   ONLY: fb_atomic_halo_list_obj, &
                                       fb_atomic_halo_obj, &
                                       fb_atomic_halo_list_get, &
                                       fb_atomic_halo_get, &
                                       fb_atomic_halo_atom_global2halo
  USE qs_fb_com_tasks_types,     ONLY: TASK_N_RECORDS, &
                                       TASK_DEST, &
                                       TASK_SRC, &
                                       TASK_PAIR, &
                                       TASK_COST, &
                                       fb_com_tasks_obj, &
                                       fb_com_atom_pairs_obj, &
                                       fb_com_atom_pairs_nullify, &
                                       fb_com_atom_pairs_has_data, &
                                       fb_com_atom_pairs_init, &
                                       fb_com_atom_pairs_create, &
                                       fb_com_atom_pairs_get, &
                                       fb_com_atom_pairs_decode, &
                                       fb_com_atom_pairs_release, &
                                       fb_com_tasks_transpose_dest_src, &
                                       fb_com_tasks_build_atom_pairs, &
                                       fb_com_tasks_set, &
                                       fb_com_tasks_get, &
                                       fb_com_tasks_nullify, &
                                       fb_com_tasks_encode_pair, &
                                       fb_com_tasks_create, &
                                       fb_com_tasks_decode_pair, &
                                       fb_com_tasks_release
  USE qs_environment_types,      ONLY: get_qs_env, &
                                       qs_environment_type
  USE cp_para_types,             ONLY: cp_para_env_type
  USE cp_dbcsr_interface,        ONLY: cp_dbcsr_type, &
                                       cp_dbcsr_get_stored_coordinates, &
                                       cp_dbcsr_get_info, &
                                       convert_sizes_to_offsets, &
                                       cp_dbcsr_get_block_p, &
                                       array_data, &
                                       cp_dbcsr_row_block_sizes, &
                                       cp_dbcsr_col_block_sizes
  USE memory_utilities,          ONLY: reallocate
  USE util,                      ONLY: sort
  USE message_passing,           ONLY: mp_alltoall
  USE termination,               ONLY: stop_program

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_atomic_matrix_methods'

  PUBLIC :: fb_atmatrix_calc_size, &
            fb_atmatrix_construct, &
            fb_atmatrix_generate_com_pairs

CONTAINS


  ! **********************************************************************
  !> \brief Calculates the atomic matrix size from a given DBCSR matrix
  !>        and atomic halo. It also calculates the first row (col) or the
  !>        row (col) atomic blocks in the atomic matrix
  !> \param dbcsr_mat : pointer to the DBCSR matrix the atomic matrix is
  !>                    to be constructed from
  !> \param atomic_halo : the atomic halo used for defining the atomic
  !>                      matrix from the DBCSR matrix
  !> \param nrows : outputs total number of rows in the atomic matrix
  !> \param ncols : outputs total number of cols in the atomic matrix
  !> \param blk_row_start : first row in each atomic blk row in the
  !>                        atomic matrix
  !> \param blk_col_start : first col in each atomic blk col in the
  !>                        atomic matrix
  !> \param error : the data container for CP2K error logs
  ! **********************************************************************
  SUBROUTINE fb_atmatrix_calc_size(dbcsr_mat, &
                                   atomic_halo, &
                                   nrows, &
                                   ncols, &
                                   blk_row_start, &
                                   blk_col_start, &
                                   error)
    TYPE(cp_dbcsr_type), POINTER :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(in) :: atomic_halo
    INTEGER, INTENT(out) :: nrows
    INTEGER, INTENT(out) :: ncols
    INTEGER, DIMENSION(:), INTENT(out) :: blk_row_start
    INTEGER, DIMENSION(:), INTENT(out) :: blk_col_start
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_calc_size', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER, DIMENSION(:), ALLOCATABLE :: first_row, last_row, &
                                          first_col, last_col
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    INTEGER :: ii, nblkrows_total, nblkcols_total, stat, natoms_in_halo

    failure = .FALSE.

    NULLIFY(halo_atoms)

    CALL cp_dbcsr_get_info(matrix=dbcsr_mat, &
                           nblkrows_total=nblkrows_total, &
                           nblkcols_total=nblkcols_total)
    ALLOCATE(first_row(nblkrows_total), last_row(nblkrows_total), &
             first_col(nblkcols_total), last_col(nblkcols_total), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL convert_sizes_to_offsets(array_data(cp_dbcsr_row_block_sizes(dbcsr_mat)), &
                                  first_row, last_row)
    CALL convert_sizes_to_offsets(array_data(cp_dbcsr_col_block_sizes(dbcsr_mat)), &
                                  first_col, last_col)
    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms, &
                            error=error)
    blk_row_start = 0
    blk_col_start = 0
    nrows = 0
    ncols = 0
    DO ii = 1, natoms_in_halo
       blk_row_start(ii) = nrows + 1
       blk_col_start(ii) = ncols + 1
       nrows = nrows + last_row(halo_atoms(ii)) - first_row(halo_atoms(ii)) + 1
       ncols = ncols + last_col(halo_atoms(ii)) - first_col(halo_atoms(ii)) + 1
    END DO
    blk_row_start(natoms_in_halo + 1) = nrows + 1
    blk_col_start(natoms_in_halo + 1) = ncols + 1
    DEALLOCATE(first_row, last_row, first_col, last_col, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  END SUBROUTINE fb_atmatrix_calc_size


  ! **********************************************************************
  !> \brief Constructs atomic matrix for filter basis method from a given
  !>        DBCSR matrix and a set of atomic send and recv pairs
  !>        corresponding to the matrix blocks that nees to be included
  !>        in the atomic matrix. The atomic matrix upon input should be
  !>        a null pointer
  ! **********************************************************************
  SUBROUTINE fb_atmatrix_construct(dbcsr_mat, &
                                   atomic_halo, &
                                   para_env, &
                                   atomic_matrix, &
                                   blk_row_start, &
                                   blk_col_start, &
                                   error)
    TYPE(cp_dbcsr_type), POINTER :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(in) :: atomic_halo
    TYPE(cp_para_env_type), POINTER :: para_env
    REAL(KIND=dp), DIMENSION(:,:), INTENT(out) :: atomic_matrix
    INTEGER, DIMENSION(:), INTENT(in) :: blk_row_start
    INTEGER, DIMENSION(:), INTENT(in) :: blk_col_start
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_construct', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure, found
    TYPE(fb_com_atom_pairs_obj) :: atom_pairs_send, atom_pairs_recv
    INTEGER(KIND=int_8), DIMENSION(:), POINTER :: pairs_send, pairs_recv
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_sizes, send_disps, &
                                          send_pair_count, send_pair_disps, &
                                          recv_sizes, recv_disps, &
                                          recv_pair_count, recv_pair_disps, &
                                          first_row, last_row, first_col, last_col
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: send_buf, recv_buf
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: mat_block
    INTEGER :: ipair, npairs_send, npairs_recv, pe, nblkrows_total, &
               nblkcols_total, arow, acol, me, ipe, ii, jj, ind, stat, &
               send_encode, recv_encode, nhalos, iatom, jatom, &
               iatom_in_halo, jatom_in_halo, numprocs, nrows_blk, ncols_blk


!BEG:DEBUG:LT:2015/01/25
    print *, "|LT| ******************** inside fb_atmatrix_construct ********************"
!END:DEBUG:LT:2015/01/25

    failure = .FALSE.

    NULLIFY(pairs_send, pairs_recv, mat_block)
    CALL fb_com_atom_pairs_nullify(atom_pairs_send)
    CALL fb_com_atom_pairs_nullify(atom_pairs_recv)

    ! initialise atomic matrix
    atomic_matrix = 0.0_dp

    ! generate send and receiv atomic pairs
    CALL fb_com_atom_pairs_create(atom_pairs_send, error)
    CALL fb_com_atom_pairs_create(atom_pairs_recv, error)
    CALL fb_atmatrix_generate_com_pairs(dbcsr_mat, &
                                        atomic_halo, &
                                        para_env, &
                                        atom_pairs_send, &
                                        atom_pairs_recv, &
                                        error)

    ! get para_env info
    numprocs = para_env%num_pe
    ! me = para_env%mepos + 1   ! my process id, starting counting from 1

    ! allocate temporary arrays
    ALLOCATE(send_sizes(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_pair_count(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_pair_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ALLOCATE(recv_sizes(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_pair_count(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_pair_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! obtain number of blocks in the DBCSR matrix and the
    ! corresponding offsets in the data_area of the DBCSR matrix
    CALL cp_dbcsr_get_info(dbcsr_mat, &
                           nblkrows_total=nblkrows_total, &
                           nblkcols_total=nblkcols_total)
    ALLOCATE(first_row(nblkrows_total), last_row(nblkrows_total), &
             first_col(nblkcols_total), last_col(nblkcols_total), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL convert_sizes_to_offsets(array_data(cp_dbcsr_row_block_sizes(dbcsr_mat)), &
                                  first_row, last_row)
    CALL convert_sizes_to_offsets(array_data(cp_dbcsr_col_block_sizes(dbcsr_mat)), &
                                  first_col, last_col)

    ! setup send buffer sizes
    CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs_send, &
                               pairs=pairs_send, &
                               npairs=npairs_send, &
                               natoms_encode=send_encode, &
                               error=error)
    send_sizes = 0
    send_pair_count = 0
    DO ipair = 1, npairs_send
       ! decode processor and (iatom, jatom) information
       CALL fb_com_atom_pairs_decode(pairs_send(ipair), &
                                     pe, iatom, jatom, send_encode)
       pe = pe + 1 ! we need proc to count from 1
       nrows_blk = last_row(iatom) - first_row(iatom) + 1
       ncols_blk = last_col(jatom) - first_col(jatom) + 1
       send_sizes(pe) = send_sizes(pe) + nrows_blk * ncols_blk
       send_pair_count(pe) = send_pair_count(pe) + 1
    END DO
    ! calculate displacements of the data of each destibation pe in
    ! send buffer and in the list of pairs to be sent
    send_disps = 0
    send_pair_disps = 0
    DO ipe = 2, numprocs
       send_disps(ipe) = send_disps(ipe-1) + send_sizes(ipe-1)
       send_pair_disps(ipe) = send_pair_disps(ipe-1) + send_pair_count(ipe-1)
    END DO
    ! allocate send buffer
    ALLOCATE(send_buf(SUM(send_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! similarly set up recieve buffer
    CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs_recv, &
                               pairs=pairs_recv, &
                               npairs=npairs_recv, &
                               natoms_encode=recv_encode, &
                               error=error)


!BEG:DEBUG:LT:2015/01/25
    print *, "|LT| fb_atmatrix_construct: npairs_send = ", npairs_send
    print *, "|LT| fb_atmatrix_construct: npairs_recv = ", npairs_recv
!END:DEBUG:LT:2015/01/25


    recv_sizes = 0
    recv_pair_count = 0
    DO ipair = 1, npairs_recv
       ! decode processor and (iatom, jatom) information
       CALL fb_com_atom_pairs_decode(pairs_recv(ipair), &
                                     pe, iatom, jatom, recv_encode)
       pe = pe + 1 ! we need proc to count from 1
       nrows_blk = last_row(iatom) - first_row(iatom) + 1
       ncols_blk = last_col(jatom) - first_col(jatom) + 1
       recv_sizes(pe) = recv_sizes(pe) + nrows_blk * ncols_blk
       recv_pair_count(pe) = recv_pair_count(pe) + 1
    END DO
    ! calculate displacements of the data of each destibation pe in
    ! recv buffer and in the list of pairs to be sent
    recv_disps = 0
    recv_pair_disps = 0
    DO ipe = 2, numprocs
       recv_disps(ipe) = recv_disps(ipe-1) + recv_sizes(ipe-1)
       recv_pair_disps(ipe) = recv_pair_disps(ipe-1) + recv_pair_count(ipe-1)
    END DO
    ! allocate recv buffer
    ALLOCATE(recv_buf(SUM(recv_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! do packing
    DO ipe = 1, numprocs
!       IF (ipe == me) CYCLE ! no need to send to self
       ! need to reuse send_sizes as an accumulative displacement, so recalculate
       send_sizes(ipe) = 0
       DO ipair = 1, send_pair_count(ipe)
          CALL fb_com_atom_pairs_decode(pairs_send(send_pair_disps(ipe) + ipair), &
                                        pe, iatom, jatom, send_encode)
          nrows_blk = last_row(iatom) - first_row(iatom) + 1
          ncols_blk = last_col(jatom) - first_col(jatom) + 1
          CALL cp_dbcsr_get_block_p(matrix=dbcsr_mat, &
                                    row=iatom, col=jatom, block=mat_block, &
                                    found=found)
          IF (.NOT. found) THEN
!TODO:LT:2014/08/07: This is not quite right here. What if the
!                    rcut-neighbor list is larger than the CP2K
!                    default neighbour_list for the matrices? Then
!                    (iatom, jatom) pair will be included in
!                    atom_pairs_send and atom_pairs_recv, but the
!                    corresponding matrix block will not be present in
!                    the DBCSR matrix.
!TODO:LT:2014/08/26: This should be corrected in the new algorithm for
!                    generating the send_pairs, which removes the pairs
!                    that does not have a corresponding DBCSR block
             CALL stop_program(routineN, moduleN, __LINE__, "Matrix block not found")
          ELSE
             ! we have found the matrix block
             DO jj = 1, ncols_blk
                DO ii = 1, nrows_blk
                   ! column major format in blocks
                   ind = send_disps(ipe) + send_sizes(ipe) + ii + (jj-1)*nrows_blk
                   send_buf(ind) = mat_block(ii,jj)
                END DO ! ii
             END DO ! jj
             send_sizes(ipe) = send_sizes(ipe) + nrows_blk*ncols_blk
          END IF
       END DO ! ipair
    END DO ! ipe

    ! do communication
    CALL mp_alltoall(send_buf, send_sizes, send_disps, &
                     recv_buf, recv_sizes, recv_disps, &
                     para_env%group)

    ! cleanup temporary arrays no longer needed
    DEALLOCATE(send_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_pair_count, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_pair_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! do unpacking
    DO ipe = 1, numprocs
!       IF (ipe == me) CYCLE
       recv_sizes(ipe) = 0
       DO ipair = 1, recv_pair_count(ipe)
          CALL fb_com_atom_pairs_decode(pairs_recv(recv_pair_disps(ipe) + ipair), &
                                        pe, iatom, jatom, recv_encode)
          nrows_blk = last_row(iatom) - first_row(iatom) + 1
          ncols_blk = last_col(jatom) - first_col(jatom) + 1
          ! get the corresponding atom indices in halo
          ! the atoms from the recv_pairs should be in the atomic_halo, because
          ! the recv_pairs are the matrix blocks requested by the local proc for
          ! this particular atomic_halo
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               iatom, iatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               jatom, jatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          ! put block into the full conventional matrix
          DO jj = 0, ncols_blk - 1
             DO ii = 0, nrows_blk - 1
                ! column major format in blocks
                ind = recv_disps(ipe) + recv_sizes(ipe) + ii + jj*nrows_blk + 1

!BEG:DEBUG:LT:2015/01/25
                print *, "|LT| fb_atmatrix_construct: ind = ", ind
                print *, "|LT| fb_atmatrix_construct: atomic_matrix(", & 
                          blk_row_start(iatom_in_halo) + ii, ",", &
                          blk_col_start(jatom_in_halo) + jj, ") = ", recv_buf(ind)
!END:DEBUG:LT:2015/01/25

                atomic_matrix(blk_row_start(iatom_in_halo) + ii, &
                              blk_col_start(jatom_in_halo) + jj) = recv_buf(ind)

             END DO ! ii
          END DO ! jj
          recv_sizes(ipe) = recv_sizes(ipe) + nrows_blk*ncols_blk
       END DO ! ipair
    END DO ! ipe

!BEG:DEBUG:LT:2015/01/25
                print *, "|LT| fb_atmatrix_construct: after assignment, atomic_matrix = "
                DO ii = 1, SIZE(atomic_matrix,1)
                   write (*,"(9f9.4)") atomic_matrix(ii,:)
                END DO
!END:DEBUG:LT:2015/01/25



    ! cleanup rest of the temporary arrays
    DEALLOCATE(recv_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_count, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(first_row, last_row, first_col, last_col, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL fb_com_atom_pairs_release(atom_pairs_send, error)
    CALL fb_com_atom_pairs_release(atom_pairs_recv, error)

  END SUBROUTINE fb_atmatrix_construct


  SUBROUTINE fb_atmatrix_generate_com_pairs(dbcsr_mat, &
                                            atomic_halo, &
                                            para_env, &
                                            atom_pairs_send, &
                                            atom_pairs_recv, &
                                            error)
    TYPE(cp_dbcsr_type), POINTER :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(in) :: atomic_halo
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(fb_com_atom_pairs_obj), INTENT(out) :: atom_pairs_send
    TYPE(fb_com_atom_pairs_obj), INTENT(out) :: atom_pairs_recv
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_generate_com_pairs', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure, found, transpose
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    TYPE(fb_com_tasks_obj) :: com_tasks_send
    TYPE(fb_com_tasks_obj) :: com_tasks_recv
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_send
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_recv
    INTEGER :: iatom, jatom, iatom_global, jatom_global, natoms_local, &
               my_basis_set_id, dest, src, ntasks_recv, ntasks_send, &
               itask, nhalos, natoms_in_halo, nblkrows_total, &
               nencode, stat, counter
    INTEGER(KIND=int_8) :: pair
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_index
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: mat_block

    failure = .FALSE.

    NULLIFY(halo_atoms, tasks_send, tasks_recv)
    CALL fb_com_tasks_nullify(com_tasks_send)
    CALL fb_com_tasks_nullify(com_tasks_recv)

    ! initialise atom_pairs_send and atom_pairs_receive
    IF (fb_com_atom_pairs_has_data(atom_pairs_send)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_send, error)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_send, error)
    END IF
    IF (fb_com_atom_pairs_has_data(atom_pairs_recv)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_recv, error)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_recv, error)
    END IF

    ! get atomic halo information
    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms, &
                            error=error)
    
    ! get the total number of atoms, we can obtain this directly
    ! from the global block row dimension of the dbcsr matrix
    CALL cp_dbcsr_get_info(matrix=dbcsr_mat, &
                           nblkrows_total=nblkrows_total)

    ! destination proc is always the local processor
    dest = para_env%mepos  ! my MPI rank

    ! generate recv task list (tasks_recv)

    ! a recv task corresponds to the copying or transfering of a
    ! matrix block in the part of the DBCSR matrix owned by the src
    ! proc to this proc in order to construct the atomic matrix
    ! corresponding to the given atomic halo. As an upper-bound, the
    ! number of matrix blocks requred do not exceed natoms_in_halo**2
    ntasks_recv = natoms_in_halo*natoms_in_halo

!BEG:DEBUG:LT:2015/01/25
    ! print *, "|LT| fb_atmatrix_generate_com_pairs: first ntasks_recv = ", ntasks_recv
!END:DEBUG:LT:2015/01/25

    ALLOCATE(tasks_recv(TASK_N_RECORDS,ntasks_recv), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! now that tasks_recv has been allocated, generate the tasks
    itask = 1
    DO iatom = 1, natoms_in_halo
       iatom_global = halo_atoms(iatom)
       DO jatom = 1, natoms_in_halo
          jatom_global = halo_atoms(jatom)
          ! find the source proc that supposed to own the block
          ! (iatom_global, jatom_global)
          transpose = .FALSE.
          CALL cp_dbcsr_get_stored_coordinates(dbcsr_mat, &
                                               iatom_global, &
                                               jatom_global, &
                                               transpose=transpose, &
                                               processor=src)
          ! we must encode the global atom indices rather the halo atomic indices
          ! in each task, because halo atomic indices are local to each halo, and
          ! each processor is working on a different halo local to them. So
          ! one processor would not have the information about the halo on another
          ! processor, rendering the halo atomic indices rather useless outside
          ! the local processor.
          tasks_recv(TASK_DEST,itask) = dest
          tasks_recv(TASK_SRC,itask) = src
          
          CALL fb_com_tasks_encode_pair(tasks_recv(TASK_PAIR,itask), &
                                        iatom_global, jatom_global, &
                                        nblkrows_total)
          ! calculation of cost not implemented at the moment
          tasks_recv(TASK_COST,itask) = 0
          itask = itask + 1
       END DO ! jatom
    END DO ! iatom

    ! create tasks
    CALL fb_com_tasks_create(com_tasks_recv, error)
    CALL fb_com_tasks_create(com_tasks_send, error)

    CALL fb_com_tasks_set(com_tasks=com_tasks_recv, &
                          task_dim=TASK_N_RECORDS, &
                          ntasks=ntasks_recv, &
                          nencode=nblkrows_total, &
                          tasks=tasks_recv, &
                          error=error)

!BEG:DEBUG:LT:2015/01/25
    ! print *, "|LT| fb_atmatrix_generate_com_pairs: ntasks_recv = ", ntasks_recv
!END:DEBUG:LT:2015/01/25

    ! genearte the send task list (tasks_send) from the recv task list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, ">", com_tasks_send, &
                                         para_env, error)


    CALL fb_com_tasks_get(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send, &
                          tasks=tasks_send, &
                          nencode=nencode, &
                          error=error)

!BEG:DEBUG:LT:2015/01/25
    ! print *, "|LT| fb_atmatrix_generate_com_pairs: ntasks_send = ", ntasks_send
!END:DEBUG:LT:2015/01/25

    ! because the atomic_halos and the neighbor_list_set used to
    ! generate the sparse structure of the DBCSR matrix do not
    ! necessarily have to coincide, we must check of the blocks in
    ! tasks_send (these should be local to the processor) do indeed
    ! exist in the DBCSR matrix, if not, then we need to prune these
    ! out of the task list

    counter = 0
    DO itask = 1, ntasks_send
       pair = tasks_send(TASK_PAIR,itask)
       CALL fb_com_tasks_decode_pair(pair, iatom_global, jatom_global, nencode)
       ! check if block exists in DBCSR matrix
       CALL cp_dbcsr_get_block_p(matrix=dbcsr_mat, &
                                 row=iatom_global, col=jatom_global, block=mat_block, &
                                 found=found)
       IF (found) THEN
          counter = counter + 1
          ! we can do this here, because essencially we are inspecting
          ! the send tasks one by one, and then omit ones which the
          ! block is not found in the DBCSR matrix. itask is always
          ! .GE. counter
          tasks_send(1:TASK_N_RECORDS,counter) = tasks_send(1:TASK_N_RECORDS,itask)
       END IF
    END DO
    ! the new send task list should have size counter. counter
    ! .LE. the old ntasks_send, thus the task list does not really
    ! need to be reallocated (as it is just a temporary array), and
    ! the useful data will cutoff at counter, and the rest of the
    ! array will just be garbage
    ntasks_send = counter

    CALL fb_com_tasks_set(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send, &
!                          tasks=tasks_send, &
!                          nencode=nencode, &
                          error=error)

    ! now, re-distribute the new send tasks list to other processors
    ! to build the updated recv tasks list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, "<", com_tasks_send, &
                                         para_env, error)


!BEG:DEBUG:LT:2015/01/25
    ! CALL fb_com_tasks_get(com_tasks=com_tasks_recv, &
    !                       ntasks=ntasks_recv, &
    !                       error=error)
    ! print *, "|LT| fb_atmatrix_generate_com_pairs: ntasks_recv = ", ntasks_recv
!END:DEBUG:LT:2015/01/25

    ! task lists are now complete, now construct the atom_pairs_send
    ! and atom_pairs_recv from the tasks lists
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_send, &
                                       atom_pairs=atom_pairs_send, &
                                       natoms_encode=nencode, &
                                       send_or_recv="send", &
                                       symmetric=.FALSE., &
                                       error=error)
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_recv, &
                                       atom_pairs=atom_pairs_recv, &
                                       natoms_encode=nencode, &
                                       send_or_recv="recv", &
                                       symmetric=.FALSE., &
                                       error=error)

    ! cleanup
    CALL fb_com_tasks_release(com_tasks_recv, error)
    CALL fb_com_tasks_release(com_tasks_send, error)

!BEG:DEBUG:LT:2015/01/20
    ! print *, "|LT| fb_atmatrix_generate_com_pairs END"
!END:DEBUG:LT:2015/01/20


  END SUBROUTINE fb_atmatrix_generate_com_pairs


  ! SUBROUTINE fb_tasks_transpose_dest_src(para_env, &
  !                                        tasks_dest_is_me, &
  !                                        ntasks_dest_is_me, &
  !                                        direction, &
  !                                        tasks_src_is_me, &
  !                                        ntasks_src_is_me, &
  !                                        error)
  !   TYPE(cp_para_env_type), POINTER :: para_env
  !   INTEGER, INTENT(inout) :: ntasks_dest_is_me
  !   INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_dest_is_me
  !   INTEGER, INTENT(inout) :: ntasks_src_is_me
  !   INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_src_is_me
  !   CHARACTER, INTENT(in) :: direction
  !   TYPE(cp_error_type), INTENT(inout) :: error

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_generate_send_tasks_from_recv_tasks', &
  !                                  routineP = moduleN//':'//routineN

  !   LOGICAL :: failure
  !   INTEGER :: itask, ipe, rank, rank_pos, task_dim, ii, jj, ind, &
  !              ntasks_in, ntasks_out
  !   INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_in, taks_out
  !   INTEGER, DIMENSION(:), ALLOCATABLE :: send_sizes, send_disps, &
  !                                         recv_sizes, recv_disps, &
  !                                         send_buf, recv_buf

  !   failure = .FALSE.

  !   NULLIFY(tasks_in, tasks_out)

  !   IF (direction == "<") THEN
  !      tasks_in => tasks_src_is_me
  !      ntasks_in = ntasks_src_is_me
  !      rank_pos = TASK_DEST
  !   ELSE
  !      tasks_in => tasks_dest_is_me
  !      ntasks_in = ntasks_dest_is_me
  !      rank_pos = TASK_SRC
  !   END IF

  !   ! allocate local arrays
  !   ALLOCATE(send_sizes(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(send_disps(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(send_buf(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   ALLOCATE(recv_sizes(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(recv_disps(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(recv_buf(para_env%num_pe), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   ! first count how many local recv tasks need to be sent to other
  !   ! processes, and share this information with the other processes.
  !   ! using send_buf as a temporary array for counting
  !   send_buf = 0
  !   ! looping over local task list
  !   DO itask = 1, ntasks_in
  !      rank = tasks_in(rank_pos,i) + 1
  !      send_buf(rank) = send_buf(rank) + 1
  !   END DO
  !   CALL mp_alltoall(send_buf, recv_buf, 1, para_env%group)

  !   ! now that we know how many recv tasks to send, pack the tasks,
  !   ! and send them around, so that the recv tasks are sent to the
  !   ! correct src processes, and these then are collected into the
  !   ! send tasks list on each of the src processes

  !   task_dim = SIZE(tasks_in, 1)
  !   send_sizes = 0
  !   send_disps = 0
  !   recv_sizes = 0
  !   recv_disps = 0

  !   ! work out the sizes of send and recv buffers and allocate them
  !   send_sizes(1) = send_buf(1) * task_dim
  !   recv_sizes(1) = recv_buf(1) * task_dim
  !   DO ipe = 2, para_env%num_pe
  !      send_sizes(ipe) = send_buf(ipe) * task_dim
  !      send_disps(ipe) = send_disps(ipe-1) + send_sizes(ipe-1)
  !      recv_sizes(ipe) = recv_buf(ipe) * task_dim
  !      recv_disps(ipe) = recv_disps(ipe-1) + recv_sizes(ipe-1)
  !   END DO

  !   ! reallocate send and recv buffers to the correct sizes for
  !   ! transferring the actual tasks
  !   DEALLOCATE(send_buf, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   DEALLOCATE(recv_buf, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(send_buf(SUM(send_sizes)), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ALLOCATE(recv_buf(SUM(recv_sizes)), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   ! now that the send buffer is of correct size, do packing
  !   send_buf = 0
  !   send_sizes = 0
  !   DO itask = 1, ntasks_in
  !      rank = tasks_in(rank_pos,itask) + 1
  !      DO ii = 1, task_dim
  !         ind = send_disps(rank) + send_sizes(rank) + ii
  !         send_buf(ind) = tasks_in(ii,itask)
  !      END DO
  !      send_sizes(rank) = send_sizes(rank) + task_dim
  !   END DO
  !   ! do communication
  !   CALL mp_alltoall(send_buf, send_sizes, send_disps, &
  !                    recv_buf, recv_sizes, recv_disps, &
  !                    para_env%group)

  !   ! deallocate send buffers
  !   DEALLOCATE(send_buf, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   DEALLOCATE(send_sizes, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   DEALLOCATE(send_disps, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   ! allocate the output task list
  !   ntasks_out = SUM(recv_sizes) / task_dim
  !   IF (direction == "<") THEN
  !      IF (ASSOCIATED(tasks_dest_is_me)) THEN
  !         DEALLOCATE(tasks_dest_is_me, STAT=stat)
  !         CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !      END IF
  !      ALLOCATE(tasks_dest_is_me(task_dim,ntasks_out), STAT=stat)
  !      tasks_out => tasks_dest_is_me
  !      ntasks_dest_is_me = ntasks_out
  !   ELSE
  !      IF (ASSOCIATED(tasks_src_is_me)) THEN
  !         DEALLOCATE(tasks_src_is_me, STAT=stat)
  !         CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !      END IF
  !      ALLOCATE(tasks_src_is_me(task_dim,ntasks_out), STAT=stat)
  !      tasks_out => tasks_src_is_me
  !      ntasks_src_is_me = ntasks_out
  !   END IF

  !   ! do unpacking
  !   itask = 0
  !   DO ipe = 1, para_env%num_pe
  !      DO ii = 0, recv_sizes(ipe) / task_dim - 1
  !         itask = itask + 1
  !         DO jj = 1, task_dim
  !            ind = recv_disps(ipe) + ii * task_dim + jj
  !            tasks_out(jj,itask) = recv_buf(ind)
  !         END DO
  !      END DO
  !   END DO

  !   ! deallocate recv buffers
  !   DEALLOCATE(recv_buf, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   DEALLOCATE(recv_sizes, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   DEALLOCATE(recv_disps, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  ! END SUBROUTINE fb_tasks_transpose_dest_src


  ! ! **********************************************************************
  ! !> \brief build atom_pair indices for the blocks from the list of tasks
  ! !> \param atom_pairs : the atom pair indices to be calclated from the
  ! !>                     tasks
  ! !> \param tasks : the list of tasks used to generate the pair indices
  ! !> \param ntasks : the total number of tasks
  ! !> \param symmetric : sets whether the matrix the blocks corresponding
  ! !>                    to the atom pairs is expected to be symmetric
  ! !> \param natoms : number of atoms corresponding to the iatom and jatom
  ! !>                 indices, this is used for decode and encode atom
  ! !>                 pair indices
  ! !> \param direction : sets if it is for send or receiv atom_pairs
  ! !> \param error : CP2K data container for error handling
  ! ! **********************************************************************
  ! SUBROUTINE fb_build_atom_pairs_from_tasks(atom_pairs, &
  !                                           tasks, &
  !                                           ntasks, &
  !                                           symmetric, &
  !                                           natoms, &
  !                                           direction, &
  !                                           error)
  !   INTEGER(KIND=int_8), DIMENSION(:), POINTER :: atom_pairs
  !   INTEGER(KIND=int_8), DIMENSION(:,:), INTENT(in) :: tasks
  !   INTEGER, INTENT(in) :: ntasks, natoms
  !   LOGICAL, INTENT(in) :: symmetric
  !   CHARACTER(LEN=*) :: direction
  !   TYPE(cp_error_type), INTENT(inout) :: error

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_build_atom_pairs_from_tasks', &
  !                                  routineP = moduleN//':'//routineN

  !   INTEGER :: itask, counter, arow, acol, pair, rank, rank_pos, stat
  !   INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_index

  !   IF (TRIM(direction) == "send") THEN
  !      rank_pos = TASK_DEST
  !   ELSE
  !      rank_pos = TASK_SRC
  !   END IF

  !   ! if atom_pairs is already assciated, reallocate to appropriate size
  !   IF (ASSOCIATED(atom_pairs)) THEN
  !      DEALLOCATE(atom_pairs, STAT=stat)
  !      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   END IF
  !   ALLOCATE(atom_pairs(ntasks), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   DO itask = 1, ntasks
  !      pair = tasks(TASK_PAIR,itask)
  !      CALL fb_decode_pair(pair, iatom, jatom, natoms)
  !      IF (symmetric) THEN
  !         IF (iatom .LE. jatom) THEN
  !            arow = iatom
  !            acol = jatom
  !         ELSE
  !            arow = jatom
  !            acol = iatom
  !         END IF
  !      ELSE
  !         arow = iatom
  !         acol = jatom
  !      END IF
  !      rank = tasks(rank_pos,itask)
  !      CALL fb_atom_pair_encode(atom_pairs(itask), &
  !                               rank, arow, acol, natoms)
  !   END DO

  !   ! sort atom_pairs so that the pairs are ordered process blocks and
  !   ! that possible duplicates may be found (we don't want to send or
  !   ! receive same information to the same destination or source more
  !   ! than once)
  !   ALLOCATE(tmp_index(ntasks), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ! only sort the actual pairs recorded in the send list
  !   CALL sort(atom_pairs, ntasks, tmp_index)
  !   DEALLOCATE(tmp_index, STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   ! remove duplicates
  !   IF (SIZE(atom_pairs) > 1) THEN
  !      counter = 1
  !      ! first atom pair must be allowed
  !      DO ii = 2, ntasks
  !         IF (atom_pairs(ii) > atom_pairs(ii-1)) THEN
  !            counter = counter + 1
  !            atom_pairs(counter) = atom_pairs(ii)
  !         END IF
  !      END DO
  !      CALL reallocate(atom_pairs, 1, counter)
  !   END IF

  ! END SUBROUTINE fb_build_atom_pairs


  ! ! **********************************************************************
  ! !> \brief encode from (iatom, jatom) pair to a single integer index
  ! !> \param int : the encoded single index
  ! !> \param iatom : iatom index of the (iatom, jatom) pair
  ! !> \param jatom : jatom index of the (iatom, jatom) pair
  ! !> \param natoms : total number of atoms corresponding to the iatom and
  ! !>                 jatom indices. This is used for encoding
  ! !> \author : Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! ! **********************************************************************
  ! SUBROUTINE fb_encode_pair(int, iatom, jatom, natoms)
  !   INTEGER(KIND=int_8), INTENT(out) :: int
  !   INTEGER, INTENT(in) :: iatom, jatom, natoms

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_encode_pair', &
  !                                  routineP = moduleN//':'//routineN

  !   INTEGER(KIND=int_8) :: natoms8

  !   int = (iatom - 1) * natoms8 + (jatom - 1)
  ! END SUBROUTINE fb_encode_pair


  ! ! **********************************************************************
  ! !> \brief decode single integer pair index to (iatom, jatom) pair
  ! !> \param int : the encoded single index
  ! !> \param iatom : iatom index of the (iatom, jatom) pair
  ! !> \param jatom : jatom index of the (iatom, jatom) pair
  ! !> \param natoms : total number of atoms corresponding to the iatom and
  ! !>                 jatom indices. This is used for encoding
  ! !> \author : Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! ! **********************************************************************
  ! SUBROUTINE fb_decode_pair(int, iatom, jatom, natoms)
  !   INTEGER(KIND=int_8), INTENT(in) :: int
  !   INTEGER, INTENT(in) :: natoms
  !   INTEGER, INTENT(out) :: iatom, jatom

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_decode_pair', &
  !                                  routineP = moduleN//':'//routineN

  !   INTEGER(KIND=int_8) :: natoms8, iatom8, jatom8

  !   natoms8 = natoms
  !   iatom8 = int / natoms8 + 1_int_8
  !   jatom8 = MOD(int, natoms8) + 1_int_8
  !   iatom = iatom8
  !   jatom = jatom8
  ! END SUBROUTINE fb_decode_pair


  ! ! **********************************************************************
  ! !> \brief encode (pe, iatom, jatom) tripplet to a single integer index
  ! !> \param ind: encoded atomic_pair_index
  ! !> \param pe: MPI process rank (to send to or recieve from) attached to
  ! !>            the pair. Note it starts counting from 0
  ! !> \param iatom: iatom index of the (iatom, jatom) pair
  ! !> \param jatom: jatom index of the (iatom, jatom) pair
  ! !> \param natoms: number of atoms corresponding to the iatom and jatom
  ! !>                indices used for encoding (pe, iatom, jatom)
  ! !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! ! **********************************************************************
  ! SUBROUTINE fb_atom_pair_encode(ind, pe, iatom, jatom, natoms)
  !   INTEGER(KIND=int_8), INTENT(out) :: ind
  !   INTEGER, INTENT(in) :: pe, iatom, jatom, natoms

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atom_pair_encode', &
  !                                  routineP = moduleN//':'//routineN

  !   INTEGER(KIND=int_8) :: natoms8, pair

  !   ! pe must start count from 0 (i.e same as MPI convension)
  !   natoms8 = natoms
  !   CALL fb_encode_pair(pair, iatom, jatom, natoms)
  !   ind = pe * natoms8 * natoms8 + pair
  ! END SUBROUTINE fb_atom_pair_encode


  ! ! **********************************************************************
  ! !> \brief decode an atom_pair index to (pe, iatom, jatom) tripplet
  ! !> \param ind: atomic_pair_index to decode
  ! !> \param pe: MPI process rank (to send to or recieve from) attached to
  ! !>            the pair. Note, starts counting from 0
  ! !> \param iatom: iatom index of the (iatom, jatom) pair
  ! !> \param jatom: jatom index of the (iatom, jatom) pair
  ! !> \param natoms: number of atoms corresponding to the iatom and jatom
  ! !>                indices used for encoding (pe, iatom, jatom) to ind
  ! !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! ! **********************************************************************
  ! SUBROUTINE fb_atom_pair_decode(ind, pe, iatom, jatom, natoms)
  !   INTEGER(KIND=int_8), INTENT(in) :: ind
  !   INTEGER, INTENT(in) :: natoms
  !   INTEGER, INTENT(out) :: pe, iatom, jatom

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atom_pair_decode', &
  !                                  routineP = moduleN//':'//routineN

  !   INTEGER(KIND=int_8) :: natoms8, pair

  !   ! pe start count from 0 (i.e same as MPI convension)
  !   natoms8 = natoms
  !   pe = ind / (natoms8*natoms8)
  !   pair =  MOD(res, natoms8*natoms8)
  !   CALL fb_decode_pair(pair, iatom, jatom, natoms)
  ! END SUBROUTINE fb_atom_pair_decode

END MODULE qs_fb_atomic_matrix_methods
