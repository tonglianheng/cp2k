MODULE qs_fb_filter_matrix_methods

  USE kinds, ONLY: dp, int_8
  USE cp_dbcsr_interface, ONLY: cp_dbcsr_type, &
                                cp_dbcsr_create, &
                                cp_dbcsr_finalize, &
                                cp_dbcsr_row_block_sizes, &
                                cp_dbcsr_col_block_sizes, &
                                array_i1d_obj, &
                                array_new, &
                                array_nullify, &
                                array_data, &
                                dbcsr_type_no_symmetry, &
                                dbcsr_distribution_obj
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_obj, &
                                     fb_atomic_halo_get, &
                                     fb_atomic_halo_list_obj, &
                                     fb_atomic_halo_list_get, &
                                     fb_atomic_halo_atom_global2halo
  USE qs_fb_trial_fns_types, ONLY: fb_trial_fns_obj, &
                                   fb_trial_fns_get
  USE qs_fb_atomic_matrix_methods, ONLY: fb_atmatrix_calc_size, &
                                         fb_atmatrix_construct
  USE qs_fb_com_tasks_types, ONLY: TASK_N_RECORDS, &
                                   TASK_DEST, &
                                   TASK_SRC, &
                                   TASK_PAIR, &
                                   TASK_COST, &
                                   fb_com_atom_pairs_obj, &
                                   fb_com_tasks_obj, &
                                   fb_com_atom_pairs_init, &
                                   fb_com_atom_pairs_create, &
                                   fb_com_atom_pairs_release, &
                                   fb_com_atom_pairs_has_data, &
                                   fb_com_atom_pairs_get, &
                                   fb_com_atom_pair_decode, &
                                   fb_com_tasks_nullify, &
                                   fb_com_tasks_encode_pair, &
                                   fb_com_tasks_create, &
                                   fb_com_tasks_set, &
                                   fb_com_tasks_transpose_dest_src, &
                                   fb_com_tasks_build_atom_pairs, &
                                   fb_com_tasks_release
  USE cp_para_types, ONLY: cp_para_env_type
  USE particle_types, ONLY: particle_type
  USE atomic_kind_types, ONLY: atomic_kind_type, &
                               get_atomic_kind
  USE cp_dbcsr_types, ONLY: cp_dbcsr_type
  USE cp_dbcsr_interface, ONLY: cp_dbcsr_get_info, &
                                cp_dbcsr_get_stored_coordinates, &
                                cp_dbcsr_put_block
  USE fermi_utils, ONLY: Fermi, FermiFixed
  ! USE qs_mo_types, ONLY: mo_set_type, &
  !                        get_mo_set
  USE lapack, ONLY: dsygv
  USE f77_blas_generic, ONLY: dgemm, &
                              dgemv
  USE string_utilities, ONLY: compress, &
                              uppercase

  IMPLICIT NONE

#include "cp_common_uses.h"

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_filter_matrix_methods'

  PUBLIC ::  fb_filrmat_build, &
             fb_fltrmat_add_blkcol, &
             fb_fltrmat_generate_com_pairs, &
             fb_fltrmat_build_atomic_fltrmat, &
             fb_fltrmat_fermi_dirac_mu, &
             fb_fltrmat_fermi_dirac_ne

CONTAINS

  SUBROUTINE fb_fltrmat_build(H_mat, &
                              S_mat, &
                              atomic_halos, &
                              trial_fns, &
                              para_env, &
                              particle_set, &
!                              mo_set, &
                              fermi_level, &
                              max_occ, &
                              filter_temp, &
                              name, &
                              filter_mat, &
                              error)
    TYPE(cp_dbcsr_type), POINTER :: H_mat
    TYPE(cp_dbcsr_type), POINTER :: S_mat
    TYPE(fb_atomic_halo_list_obj), INTENT(in) :: atomic_halos
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
!    TYPE(mo_set_type), POINTER :: mo_set
    REAL(KIND=dp), INTENT(in) :: fermi_level
    REAL(KIND=dp), INTENT(in) :: max_occ
    REAL(KIND=dp), INTENT(in) :: filter_temp
    CHARACTER(LEN=*), INTENT(in) :: name
    TYPE(cp_dbcsr_type), POINTER :: filter_mat
    TYPE(cp_error_type), INTENT(inout) :: error    

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_fltrmat_build', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    CHARACTER(LEN=default_string_length) :: name_string
    CHARACTER(LEN=32) :: symmetry_string
    TYPE(dbcsr_distribution_obj), POINTER :: dbcsr_dist
    TYPE(array_i1d_obj) :: col_blk_size, row_blk_size
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER, DIMENSION(:), POINTER :: col_blk_size_data, ntfns
    INTEGER :: ikind, ihalo, iblkcol, nblkcols_total, nhalos, stat
    TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
    

    failure = .FALSE.

    NULLIFY(halos, atomic_kind, col_blk_size_data, dbcsr_dist, ntfns)
    CALL array_nullify(row_blk_size)
    CALL array_nullify(col_blk_size)

    ! filter_mat must be of a dissassociated status (i.e. brand new)
    CPPrecondition(.NOT.ASSOCIATED(filter_mat), cp_failure_level, routineP, error, failure)
    
    ! get trial function information
    CALL fb_trial_fns_get(trial_fns=trial_fns, &
                          nfunctions=ntfns, &
                          error=error)

    ! calculate the row_blk_size and col_blk_size arrays for
    ! constructing the filter matrix in DBCSR format
    ! row_blk_size for the filter matrix is the same as H or S
    CALL cp_dbcsr_get_info(mat_H%matrix, &
                           nblkcols_total=nblkcols_total, &
                           row_blk_size=row_blk_size, &
                           distribution=dbcsr_dist)
    ALLOCATE(col_blk_size_data(nblkcols_total), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)    
    col_blk_size_data = 0
    DO iblkcol = 1, nblkcols_total
       atomic_kind => particle_set(iblkcol)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=ikind)
       col_blk_size_data(iblkcol) = ntfns(ikind)
    END DO
    ! DO NOT deallocate cbs if gift=.TRUE. as col_blk_sizes will only point to cbs
    CALL array_new(col_blk_size, col_blk_size_data, gift=.TRUE.)
    name_string = name
    CALL compress(name_string)
    CALL uppercase(name_string)
    ! the filter matrix is non-square and is always non-symmetric
    symmetry_string = dbcsr_type_no_symmetry
    ! create empty filter matrix
    ALLOCATE(filter_mat, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL cp_dbcsr_init(filter_mat, error=error)
    CALL cp_dbcsr_create(matrix=filter_mat, &
                         name=name_string, &
                         dist=dbcsr_dist, &
                         matrix_type=symmetry_string, &
                         row_blk_size=row_blk_size, &
                         col_blk_size=col_blk_size, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
                                 nhalos=nhalos, &
                                 halos=halos, &
                                 error=error)
    DO ihalo = 1, nhalos
       ! construct the filter matrix block by block
       CALL fb_fltrmat_add_blkcol(H_mat, &
                                  S_mat, &
                                  halos(ihalo), &
                                  trial_fns, &
                                  para_env, &
!                                  mo_set, &
                                  fermi_level, &
                                  max_occ, &
                                  filter_temp, &
                                  filter_mat, &
                                  error)
    END DO
    ! finalise the filter matrix
    CALL cp_dbcsr_finalize(filter_mat, error=error)

  END SUBROUTINE fb_fltrmat_build
  
  

  SUBROUTINE fb_fltrmat_add_blkcol(H_mat, &
                                   S_mat, &
                                   atomic_halo, &
                                   trial_fns, &
                                   para_env, &
                                   particle_set, &
                                   fermi_level, &
                                   max_occ, &
!                                   mo_set, &
                                   filter_temp, &
                                   filter_mat, &
                                   error)
    TYPE(cp_dbcsr_type), POINTER :: H_mat
    TYPE(cp_dbcsr_type), POINTER :: S_mat
    TYPE(fb_atomic_halo_obj), INTENT(in) :: atomic_halo
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
!    TYPE(mo_set_type), POINTER :: mo_set
    REAL(KIND=dp), INTENT(in) :: fermi_level
    REAL(KIND=dp), INTENT(in) :: max_occ
    REAL(KIND=dp), INTENT(in) :: filter_temp
    TYPE(cp_dbcsr_type), POINTER :: filter_mat
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_fltrmat_add_blkcol', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(fb_com_atom_pairs_obj) :: com_pairs_send, com_pairs_recv
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER :: ipair, stat, pe, ipe, numprocs &
               iatom_global, jatom_global, &
               iatom_in_halo, jatom_in_halo, &
               nrows_blk, ncols_blk, &
               npairs_send, send_encode, send_sizes, send_pair_count, &
               npairs_recv, recv_encode, recv_sizes, recv_pair_count, &
               natoms_in_halo, nrows_atmatrix, ncols_atmatrix, &
               mat_dim, jkind, itrial, ind
    INTEGER(KIND=int_8), DIMENSION(:), POINTER :: pairs_send, pairs_recv
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_sizes, send_disps, &
                                          send_pair_count, send_pair_disps, &
                                          recv_sizes, recv_disps, &
                                          recv_pair_count, recv_pair_disps, &
                                          first_row, last_row, first_col, last_col, &
                                          atomic_H_blk_row_start, &
                                          atomic_H_blk_col_start, &
                                          atomic_S_blk_row_start, &
                                          atomic_S_blk_col_start
    INTEGER, DIMENSION(:), POINTER :: ntfns
    INTEGER, DIMENSION(:,:), POINTER :: tfns
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: send_buf, recv_buf, &
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: atomic_H, atomic_S
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE ::  atomic_filter_mat

    failure = .FALSE.

    NULLIFY(pairs_send, pairs_recv, ntfns, tfns, atomic_kind)
    CALL fb_com_atom_pairs_nullify(com_pairs_send, error)
    CALL fb_com_atom_pairs_nullify(com_pairs_recv, error)


    ! ----------------------------------------------------------------------
    ! Get communication buffers ready
    ! ----------------------------------------------------------------------

    ! generate send and recv atom pairs
    CALL fb_com_atom_pairs_create(atom_pairs_send, error)
    CALL fb_com_atom_pairs_create(atom_pairs_recv, error)
    CALL fb_fltrmat_generate_com_pairs(filter_mat, &
                                       atomic_halo, &
                                       para_env, &
                                       com_pairs_send, &
                                       com_pairs_recv, &
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

    ! obtain information from the DBCSR filter matrix
    CALL cp_dbcsr_get_info(filter_mat, &
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

    ! obtain trail function information
    CALL fb_trial_fns_get(trial_fns=trial_fns, &
                          nfunctions=ntfns, &
                          functions=tfns, &
                          error)

    ! setup send buffer sizes
    CALL fb_com_atom_pairs_get(atom_pairs=com_pairs_send, &
                               pairs=pairs_send, &
                               npairs=npairs_send, &
                               natoms_encode=send_encode, &
                               error=error)
    send_sizes = 0
    send_pair_count = 0
    DO ipair = 1, npairs_send
       ! decode processor and (iatom, jatom) information
       CALL fb_com_atom_pair_decode(pairs_send(ipair), &
                                    pe, iatom_global, jatom_global, &
                                    send_encode)
       atomic_kind => particle_set(jatom_global)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=jkind)
       nrows_blk = last_row(iatom_global) - first_row(iatom_global) + 1
       ncols_blk = ntfns(jkind)
       pe = pe + 1 ! we need proc to count from 1
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
    CALL fb_com_atom_pairs_get(atom_pairs=com_pairs_recv, &
                               pairs=pairs_recv, &
                               npairs=npairs_recv, &
                               natoms_encode=recv_encode, &
                               error=error)
    recv_sizes = 0
    recv_pair_count = 0
    DO ipair = 1, npairs_recv
       ! decode processor and (iatom, jatom) information
       CALL fb_com_atom_pair_decode(pairs_recv(ipair), &
                                    pe, iatom_global, jatom_global, &
                                    recv_encode)
       atomic_kind => particle_set(jatom_global)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=jkind)
       nrows_blk = last_row(iatom_global) - first_row(iatom_global) + 1
       ncols_blk = ntfns(jkind)
       pe = pe + 1 ! we need proc to count from 1
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

    ! ----------------------------------------------------------------------
    ! Construct atomic filter matrix for this atomic_halo
    ! ----------------------------------------------------------------------

    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            error=error)

    ! construct atomic matrix for H for atomic_halo
    ALLOCATE(atomic_H_blk_row_start(natoms_in_halo + 1), &
             atomic_H_blk_col_start(natoms_in_halo + 1), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_calc_size(mat_H, &
                               atomic_halo, &
                               nrows_atmatrix, &
                               ncols_atmatrix, &
                               atomic_H_blk_row_start, &
                               atomic_H_blk_col_start, &
                               error)
    ALLOCATE(atomic_H(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_construct(mat_H, &
                               atomic_halo, &
                               para_env, &
                               atomic_H, &
                               atomic_H_blk_row_start, &
                               atomic_H_blk_col_start, &
                               error)

    ! construct atomic matrix for S for atomic_halo
    ALLOCATE(atomic_S_blk_row_start(natoms_in_halo + 1), &
             atomic_S_blk_col_start(natoms_in_halo + 1), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_calc_size(mat_S, &
                               atomic_halo, &
                               nrows_atmatrix, &
                               ncols_atmatrix, &
                               atomic_S_blk_row_start, &
                               atomic_S_blk_col_start, &
                               error)
    ALLOCATE(atomic_S(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_construct(mat_S, &
                               atomic_halo, &
                               para_env, &
                               atomic_S, &
                               atomic_S_blk_row_start, &
                               atomic_S_blk_col_start, &
                               error)

    ! construct the atomic filter matrix
    ALLOCATE(atomic_filter_mat(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_fltrmat_build_atomic_fltrmat(atomic_H, &
                                         atomic_S, &
                                         fermi_level, &
                                         max_occ, &
                                         filter_temp, &
                                         atomic_filter_mat, &
                                         error)

    ! ----------------------------------------------------------------------
    ! Construct filter matrix blocks and add to the correct locations
    ! in send_buffer
    ! ----------------------------------------------------------------------

    ! initialise send buffer
    send_buf = 0.0_dp
    ! assign values
    DO ipe = 1, numprocs
       send_sizes(ipe) = 0
       DO ipair = 1, send_pair_count(ipe)
          CALL fb_com_atom_pair_decode(pairs_send(send_pair_disps(ipe) + ipair), &
                                       pe, iatom_global, jatom_global, &
                                       send_encode)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               iatom_global, iatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               jatom_global, jatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          atomic_kind => particle_set(jatom_global)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, &
                               kind_number=jkind)
          nrows_blk = last_row(iatom_global) - first_row(iatom_global) + 1
          ncols_blk = ntfns(jkind)
          ! do it column-wise one trial function at a time
          DO itrial = 1, ntfns(jkind)
             ind = send_disp(ipe) + send_sizes(ipe) + (itrial-1) * nrows_blk
             CALL dgemv(TRANSA = "N",                                          &
                        TRANSB = "N",                                          &
                        M      = nrows_blk,                                    &
                        N      = ncols_atmatrix,                               &
                        ALPHA  = 1.0_dp,                                       &
                        A      = atomic_filter_mat(                            &
                                   atomic_H_blk_row_start(iatom_in_halo):, :   &
                                   ! atomic_H_blk_row_start(iatom_in_halo + 1) - 1, &
                                   ! 1 : ncols_atmatrix                          &
                                 ),                                            &
                        LDA    = nrows_atmatrix,                               &
                        X      = atomic_S(                                     &
                                   1 : nrows_atmatrix,                         &
                                   atomic_S_blk_row_start(jhalo_atom) +        &
                                   tfns(itrial,jkind)                          &
                                 ),                                            &
                        INCX   = 1,                                            &
                        BETA   = 0.0_dp,                                       &
                        Y      = send_buf(ind + 1 : ind + nrows_blk),          &
                        INCY   = 1)
          END DO ! itrial
          send_sizes(ipe) = send_sizes(ipe) + nrows_blk * ncols_blk
       END DO ! ipair
    END DO  ! ipe

    DEALLOCATE(atomic_H, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_H_blk_row_start, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_S, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_S_blk_row_start, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_filter_mat, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! ----------------------------------------------------------------------
    ! Do communication
    ! ----------------------------------------------------------------------

    CALL mp_alltoall(send_buf, send_sizes, send_disps, &
                     recv_buf, recv_sizes, recv_disps, &
                     para_env&group)

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
    DEALLOCATE(recv_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)


    ! ----------------------------------------------------------------------
    ! Unpack the recv buffer and add the blocks to correct parts of
    ! the DBCSR filter matrix
    ! ----------------------------------------------------------------------

    DO ipe = 1, numprocs
       recv_sizes(ipe) = 0
       DO ipair = 1, recv_pair_count(ipe)
          CALL fb_com_atom_pair_decode(pairs_recv(recv_pair_disps(ipe) + ipair), &
                                       pe, iatom_global, jatom_global, &
                                       send_encode)
          atomic_kind => particle_set(jatom_global)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, &
                               kind_number=jkind)
          nrows_blk = last_row(iatom_global) - first_row(iatom_global) + 1
          ncols_blk = ntfns(jkind)
          ind = recv_disp(ipe) + recv_sizes(ipe)
          CALL cp_dbcsr_put_block(filter_mat, &
                                  iatom_global, jatom_global, &
                                  recv_buff(ind + 1 : ind + nrows_blk * ncols_blk))
       END DO ! ipair
       recv_sizes(ipe) = recv_sizes(ipe) + nrows_blk * ncols_blk
    END DO ! ipe

    ! cleanup rest of the temporary arrays
    DEALLOCATE(recv_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_count, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(first_row, last_row, first_col, last_col, STAT)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL fb_com_atom_pairs_release(atom_pairs_send, error)
    CALL fb_com_atom_pairs_release(atom_pairs_recv, error)

    ! cannot finalise the matrix until all blocks has been added

  END SUBROUTINE fb_fltrmat_add_blkcol



  SUBROUTINE fb_fltrmat_generate_com_pairs(filter_mat, &
                                           atomic_halo, &
                                           para_env, &
                                           atom_pairs_send, &
                                           atom_pairs_recv, &
                                           error)
    TYPE(cp_dbcsr_type), POINTER :: filter_mat
    TYPE(fb_atomic_halo_obj), INTENT(in) :: atomic_halo
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(fb_com_atom_pairs_obj), INTENT(out) :: atom_pairs_send
    TYPE(fb_com_atom_pairs_obj), INTENT(out) :: atom_pairs_recv
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_fltrmat_generate_com_tasks', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(fb_com_tasks_obj) :: com_tasks_send, com_tasks_recv
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_send, tasks_recv
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    INTEGER :: itask, src, dest, jatom_global, iatom_in_halo, &
               natoms_in_halo, ntasks_send, ntasks_recv, nblkrows_total, &
               stat

    failure = .FALSE.

    NULLIFY(tasks_send, tasks_recv)
    CALL fb_com_tasks_nullify(com_tasks_send, error)
    CALL fb_com_tasks_nullify(com_tasks_recv, error)

    ! initialise atom_pairs_send and atom_pairs_recv
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

    ! source is always the local processor
    src = para_env%mepos

    ! The total number of filter matrix blocks each processor is going
    ! to construct equals to the total number of halo atoms in all of
    ! the atomic halos local to the processor. The number of send
    ! tasks will not exceed this. We do one halo (col) at a time, and
    ! each call of this subroutine will only work on one filter matrix
    ! col corresponding to atomic_halo.

    ! The col atom block index for each filter matrix block are the
    ! owner atom of each halo. The row atom block index for each
    ! filter matrix block corresponding to each col are the halo atoms
    ! of the corresponding halos. Filter matrix is non-symmetric: it
    ! is non-square, because the blocks themselves are non-sqaure

    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            owner_atom=jatom_global, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms, &
                            error=error)
    ntasks_send = natoms_in_halo

    ! allocate send tasks
    ALLOCATE(tasks_send(TASK_N_RECORDS,ntasks_send), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! Get the total number of atoms, this can be obtained from the
    ! total number of block rows in the DBCSR filter matrix.  We
    ! assumes that before calling this subroutine, the filter_mat has
    ! already been created and initialised: i.e. using
    ! cp_dbcsr_create_new. Even if the matrix is at the moment empty,
    ! the attribute nblkrows_total is already assigned from the dbcsr
    ! distribution data
    CALL cp_dbcsr_get_info(filter_mat, &
                           nblkrows_total=nblkrows_total)

    ! construct send tasks
    itask = 1
    DO iatom_in_halo = 1, natoms_in_halo
       iatom_global = halo_atoms(iatom_in_halo)
       iatom_stored = iatom_global
       jatom_stored = jatom_global
       ! find where the constructed block of filter matrix belongs to
       CALL cp_dbcsr_get_stored_coordinates(filter_mat, &
                                            iatom_stored, &
                                            jatom_stored, &
                                            transpose=.FALSE., &
                                            processor=dest)
       ! create the send tasks
       tasks_send(TASK_DEST,itask) = dest
       tasks_send(TASK_SRC,itask) = src
       CALL fb_com_tasks_encode_pair(tasks_send(TASK_PAIR,itask), &
                                     iatom_global, jatom_global, &
                                     nblkrows_total)
       ! calculation of cost not implemented at the moment
       tasks_send(TASK_COST,itask) = 0
       itask = itask + 1
    END DO ! iatom_in_halo

    ! construct send tasks
    CALL fb_com_tasks_create(com_tasks_recv, error)
    CALL fb_com_tasks_create(com_tasks_send, error)

    CALL fb_com_tasks_set(com_tasks=com_tasks_send, &
                          task_dim=TASK_N_RECORDS, &
                          ntasks=ntasks_send, &
                          nencode=nblkrows_total, &
                          tasks=tasks_send, &
                          error=error)

    ! genearte the recv task list (tasks_recv) from the send task list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, "<", com_tasks_send, &
                                         para_env, error)

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

  END SUBROUTINE fb_fltrmat_generate_com_pairs


  SUBROUTINE fb_fltrmat_build_atomic_fltrmat(atomic_H, &
                                             atomic_S, &
!                                             mo_set, &
                                             fermi_level, &
                                             max_occ, &
                                             filter_temp, &
                                             atomic_filter_mat, &
                                             error)
    REAL(KIND=dp), DIMENSION(:,:), INTENT(inout) :: atomic_H
    REAL(KIND=dp), DIMENSION(:,:), INTENT(in) :: atomic_S
 !   TYPE(mo_set_type), POINTER :: mo_set
    REAL(KIND=dp), INTENT(in) :: fermi_level
    REAL(KIND=dp), INTENT(in) :: max_occ
    REAL(KIND=dp), INTENT(in) :: filter_temp
    REAL(KIND=dp), DIMENSION(:,:), INTENT(out) :: atomic_filter_mat
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_fltrmat_build_atomic_fltrmat', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat, info, mat_dim, work_array_size

!    REAL(KIND=dp) :: fermi_level, maxocc
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: filter_function, &
                                                eigenvalues, work
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: filtered_eigenvectors


    failure = .FALSE.

    ! initialise
    atomic_filter_mat = 0.0_dp
    mat_dim = SIZE(atomic_H, 1)

    ! diagonalise using LAPACK
    ALLOCATE(eigenvalues(mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! get optimal work array size
    ALLOCATE(work(1), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL dsygv(ITYPE = 1,           &
               JOBZ  = 'V',         &
               UPLO  = 'U',         &
               N     = mat_dim,     &
               A     = atomic_H,    &
               LDA   = mat_dim,     &
               B     = atomic_S,    &
               LDB   = mat_dim,     &
               W     = eigenvalues, &
               WORK  = work,        &
               LWORK = -1,          &
               INFO  = info)
    work_array_size = work(1)
    ! now allocate work array
    DEALLOCATE(work, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(work(work_array_size), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! do calculation
    CALL dsygv(ITYPE = 1,               &
               JOBZ  = 'V',             &
               UPLO  = 'U',             &
               N     = mat_dim,         &
               A     = atomic_H,        &
               LDA   = mat_dim,         &
               B     = atomic_S,        &
               LDB   = mat_dim,         &
               W     = eigenvalues,     &
               WORK  = work,            &
               LWORK = work_array_size, &
               INFO  = info)
    ! check if diagonalisation is successful
    IF (info .NE. 0) THEN
       WRITE (*,*) "DSYGV ERROR MESSAGE: ", info
       CPErrorMessage(cp_failure_level, routineP, "DSYGV failed", error)
       CPPrecondition(.FALSE., cp_failure_level, routineP, error, failure)
    END IF

    DEALLOCATE(work, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! atomic_H now contains the eigenvectors

    ! first get the filter function
    ALLOCATE(filter_function(mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! ! get system mu
    ! CALL get_mo_set(mo_set=mo_set, &
    !                 mu=fermi_level, &
    !                 maxocc=maxocc)
    CALL fb_fltrmat_fermi_dirac_mu(filter_function, &
                                   eigenvalues, &
                                   filter_temp, &
                                   fermi_level, &
                                   max_occ)

    DEALLOCATE(eigenvalues, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! atomic_H has the eigenvectors, construct the version of it
    ! filtered through the filter function
    ALLOCATE(filtered_eigenvectors(mat_dim,mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO jj = 1, mat_dim
       DO ii = 1, mat_dim
          filtered_eigenvectors(ii,jj) = &
               filter_function(jj) * atomic_H(ii,jj)
       END DO ! ii
    END DO ! jj

    DEALLOCATE(filter_function, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! construct atomic filter matrix
    CALL dgemm(TRANSA = "N",                   &
               TRANSB = "T",                   &
               M      = mat_dim,               &
               N      = mat_dim,               &
               K      = mat_dim,               &
               ALPHA  = 1.0_dp,                &
               A      = filtered_eigenvectors, &
               LDA    = mat_dim,               &
               B      = atomic_H,              &
               LDB    = mat_dim,               &
               BETA   = 0.0_dp,                &
               C      = atomic_filter_mat,     &
               LDC    = mat_dim)

    DEALLOCATE(filtered_eigenvectors, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE fb_fltrmat_build_atomic_fltrmat


  SUBROUTINE fb_fltrmat_fermi_dirac_mu(f, eigenvals, T, mu, maxocc)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in) :: T
    REAL(KIND=dp), INTENT(in) :: mu
    REAL(KIND=dp), INTENT(in) :: maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_fltrmat_fermi_dirac_mu', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: ne
    REAL(KIND=dp) :: kTS

    CALL Fermi(f, ne, kTS, eigenvals, mu, T, maxocc)
  END SUBROUTINE fb_fltrmat_fermi_dirac_mu


  SUBROUTINE fb_fltrmat_fermi_dirac_ne(f, eigenvals, T, ne, maxocc)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in) :: T
    REAL(KIND=dp), INTENT(in) :: ne
    REAL(KIND=dp), INTENT(in) :: maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_fltrmat_fermi_dirac_ne', &
                                   routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: mu, kTS
    ! mu is the calculated fermi level
    ! kTS is the calculated entropic contribution to the energy i.e. -TS
    ! kTS= kT*[f ln f + (1-f) ln (1-f)]
    LOGICAL       :: failure

    failure = .FALSE.
    CALL FermiFixed(f, mu, kTS, eigenval, ne, fermi_T, maxocc)
  END SUBROUTINE fb_fltrmat_fermi_dirac_ne


END MODULE qs_fb_filter_matrix_methods
