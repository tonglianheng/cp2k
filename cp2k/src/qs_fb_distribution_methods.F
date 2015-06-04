!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

Module qs_fb_distribution_methods

  USE kinds, ONLY: dp
  USE cell_types, ONLY: cell_type
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE particle_types, ONLY: particle_type
  USE qs_fb_atomic_halo_types, ONLY: &
       fb_atomic_halo_obj, fb_atomic_halo_set, fb_build_pair_radii, &
       fb_atomic_halo_create, fb_atomic_halo_init, &
       fb_atomic_halo_build_halo_atoms, fb_atomic_halo_cost, &
       fb_atomic_halo_release, fb_atomic_halo_nullify
  USE qs_fb_env_types, ONLY: fb_env_obj,&
                             fb_env_get,&
                             fb_env_set
  USE qs_kind_types, ONLY: qs_kind_type
  USE cp_dbcsr_interface, ONLY: cp_dbcsr_get_info,&
                                cp_dbcsr_type,&
                                cp_dbcsr_p_type,&
                                cp_dbcsr_nblkcols_total,&
                                dbcsr_distribution_obj,&
                                dbcsr_distribution_mp,&
                                dbcsr_mp_obj,&
                                dbcsr_mp_nprows,&
                                dbcsr_mp_npcols,&
                                dbcsr_mp_get_process
  USE cp_para_types, ONLY: cp_para_env_type
  USE util, ONLY: sort
  USE timings, ONLY: timeset,&
                     timestop
  USE input_section_types, ONLY: section_vals_type
  USE cp_output_handling, ONLY: cp_print_key_finished_output,&
                                cp_print_key_unit_nr                                            

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fb_distribution_build

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_distribution_methods'

  TYPE fb_distribution_element
     INTEGER :: id
     REAL(KIND=dp) :: cost
  END type fb_distribution_element

  TYPE fb_distribution_list
     TYPE(fb_distribution_element), DIMENSION(:), POINTER :: list => NULL()
     INTEGER :: nelements
     REAL(KIND=dp) :: cost
  END type fb_distribution_list

  TYPE fb_preferred_procs_list
     INTEGER, DIMENSION(:), POINTER :: list => NULL()
     INTEGER :: nprocs
  END type fb_preferred_procs_list

! Parameters related to automatic resizing of the hash_table:
! Resize by EXPAND_FACTOR if total no. slots / no. of filled slots < ENLARGE_RATIO
  INTEGER, PARAMETER, PRIVATE :: ENLARGE_RATIO = 1
  INTEGER, PARAMETER, PRIVATE :: REDUCE_RATIO = 3
  INTEGER, PARAMETER, PRIVATE :: EXPAND_FACTOR = 2
  INTEGER, PARAMETER, PRIVATE :: SHRINK_FACTOR = 2

  INTERFACE fb_distribution_remove
     MODULE PROCEDURE fb_distribution_remove_ind, &
                      fb_distribution_remove_el
  END INTERFACE fb_distribution_remove

  INTERFACE fb_distribution_move
     MODULE PROCEDURE fb_distribution_move_ind, &
                      fb_distribution_move_el
  END INTERFACE fb_distribution_move

CONTAINS

  SUBROUTINE fb_distribution_build(fb_env, qs_env, scf_section, error)
    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_build', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure, move_happened, acceptable_move
    INTEGER :: natoms, nkinds, nprocs, my_pe, stat, iatom, &
         owner_id_in_halo, nhalo_atoms, n_common_sets, ipe, &
         lowest_cost_ind, ii, i_common_set, handle, lb, ub, &
         pref_pe
    REAL(KIND=dp) :: average_cost
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: mat_ks
    REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: pair_radii
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), POINTER :: qs_kind_set
    TYPE(fb_atomic_halo_obj) :: atomic_halo
    INTEGER, DIMENSION(:), POINTER :: halo_atoms, local_atoms
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: cost_per_atom, &
         cost_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: common_set_ids, &
         pos_in_preferred_list, pe, local_atoms_all, &
         local_atoms_starts, local_atoms_sizes
    TYPE(fb_preferred_procs_list), DIMENSION(:), ALLOCATABLE :: &
         preferred_procs_set
    TYPE(fb_distribution_list), DIMENSION(:), ALLOCATABLE :: dist
    TYPE(fb_distribution_element) :: element

    CALL timeset(routineN, handle)

    failure = .FALSE.
    
    NULLIFY(mat_ks, rcut, cell, para_env, particle_set, qs_kind_set, &
            halo_atoms, local_atoms)
    CALL fb_atomic_halo_nullify(atomic_halo)

    ! obtain relevant data from fb_env, qs_env
    CALL fb_env_get(fb_env=fb_env, &
                    rcut=rcut, &
                    error=error)
    CALL get_qs_env(qs_env=qs_env, &
                    natom=natoms, &
                    particle_set=particle_set, &
                    qs_kind_set=qs_kind_set, &
                    nkind=nkinds, &
                    cell=cell, &
                    para_env=para_env, &
                    matrix_ks=mat_ks, &
                    error=error)
    nprocs = para_env%num_pe
    my_pe = para_env%mepos + 1 ! counting from 1

    ! for each global atom, build atomic halo and get the associated cost
    ALLOCATE(pair_radii(nkinds,nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_build_pair_radii(rcut, nkinds, pair_radii)
    CALL fb_atomic_halo_create(atomic_halo, error)
    ALLOCATE(cost_per_atom(natoms), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO iatom = 1, natoms
       CALL fb_atomic_halo_init(atomic_halo, error)
       CALL fb_atomic_halo_build_halo_atoms(iatom, &
                                            particle_set, &
                                            cell, &
                                            pair_radii, &
                                            halo_atoms, &
                                            nhalo_atoms, &
                                            owner_id_in_halo, &
                                            error)
       CALL fb_atomic_halo_set(atomic_halo=atomic_halo, &
                               owner_atom=iatom, &
                               natoms=nhalo_atoms, &
                               halo_atoms=halo_atoms, &
                               error=error)

!BEG:DEBUG:LT:2015/06/04
       write (my_pe+200,*) "atomic halo for atom: ", iatom
       DO ii = 1, SIZE(halo_atoms)
          write (my_pe+200,*) halo_atoms(ii)
       END DO
!END:DEBUG:LT:2015/06/04

       NULLIFY(halo_atoms)
       cost_per_atom(iatom) = fb_atomic_halo_cost(atomic_halo, particle_set, qs_kind_set)
    END DO
    DEALLOCATE(pair_radii, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atomic_halo_release(atomic_halo, error)

    ! build the preferred_procs_set according to DBCSR mat H
    ALLOCATE(preferred_procs_set(natoms), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(common_set_ids(natoms), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_build_preferred_procs(mat_ks(1)%matrix, &
                                  natoms, &
                                  preferred_procs_set, &
                                  common_set_ids, &
                                  n_common_sets, &
                                  error)

    ! for each atomic halo, construct distribution_element, and assign
    ! the element to a processors using preferred_procs_set in a
    ! round-robin manner
    ALLOCATE(dist(nprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ipe = 1, nprocs
       CALL fb_distribution_init(dist=dist(ipe), error=error)
    END DO
    ALLOCATE(pos_in_preferred_list(n_common_sets), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    pos_in_preferred_list(:) = 0
    DO iatom = 1, natoms
       element%id = iatom
       element%cost = cost_per_atom(iatom)
       i_common_set = common_set_ids(iatom)
       pos_in_preferred_list(i_common_set) = &
            MOD(pos_in_preferred_list(i_common_set), &
                preferred_procs_set(iatom)%nprocs) + 1
       ipe = preferred_procs_set(iatom)%list(pos_in_preferred_list(i_common_set))

!BEG:DEBUG:LT:2015/06/03
       write (my_pe+200,*) "|LT| build: iatom = ", iatom
       write (my_pe+200,*) "|LT| build: preferred_procs_set = ", preferred_procs_set(iatom)%list(:)
       write (my_pe+200,*) "|LT| build: ipe = ", ipe
!END:DEBUG:LT:2015/06/03

       CALL fb_distribution_add(dist(ipe), element, error)

!BEG:DEBUG:LT:2015/06/03
       write (my_pe+200,*) "|LT| dist(",ipe,") = "
       DO ii = 1, dist(ipe)%nelements
          write (my_pe+200,*) dist(ipe)%list(ii)%id, dist(ipe)%list(ii)%cost
       END DO
       write (my_pe+200,*) "|LT|   cost for this dist = ", dist(ipe)%cost
!END:DEBUG:LT:2015/06/03

    END DO

!BEG:DEBUG:LT:2015/06/04
    CALL fb_distribution_write_info(dist, scf_section, error)
!END:DEBUG:LT:2015/06/04

    DEALLOCATE(pos_in_preferred_list, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(common_set_ids, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(cost_per_atom, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! sort processors according to the overall cost of their assigned
    ! corresponding distribution
    ALLOCATE(cost_per_proc(nprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ipe = 1, nprocs
       cost_per_proc(ipe) = dist(ipe)%cost
    END DO
    ALLOCATE(pe(nprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL sort(cost_per_proc, nprocs, pe)
    ! now that cost_per_proc is sorted, ipe's no longer give mpi
    ! ranks, the correct one to use should be pe(ipe)

    ! work out the ideal average cost per proc if work load is evenly
    ! distributed
    average_cost = SUM(cost_per_proc) / REAL(nprocs, dp)

    DEALLOCATE(cost_per_proc, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! loop over the processors, starting with the highest cost, move
    ! atoms one by one:
    !   1. FIRST to the next processor in the preferred list that has
    !      cost below average. IF no such proc is found, THEN
    !   2. to the next procesor in the overall list that has cost
    !      below average.
    ! repeat until the cost on this processor is less than or equal
    ! to the average cost
    lowest_cost_ind = 1
    DO ipe = nprocs, 1, -1
       redistribute: DO WHILE (dist(pe(ipe))%cost .GT. average_cost)
          iatom = dist(pe(ipe))%list(lowest_cost_ind)%id
          move_happened = .FALSE.
          ! first try to move to a preferred process
          preferred: DO ii = 1, preferred_procs_set(iatom)%nprocs
             pref_pe = preferred_procs_set(iatom)%list(ii)
             acceptable_move = &
                  fb_distribution_acceptable_move(dist(pe(ipe)), &
                                                  dist(pe(ipe))%list(lowest_cost_ind), &
                                                  dist(pref_pe), &
                                                  average_cost)
             IF ((pref_pe .NE. pe(ipe)) .AND. acceptable_move) THEN
                CALL fb_distribution_move(dist(pe(ipe)), &
                                          lowest_cost_ind, &
                                          dist(pref_pe), &
                                          error)
                move_happened = .TRUE.
                EXIT preferred
             END IF
          END DO preferred
          ! if no preferred process is avaliable, move to a proc in
          ! the sorted list that has cost less than average.  remember
          ! that some of the proc may have already taken redistributed
          ! atoms, and thus may beome unavaliable (full)
          IF (.NOT. move_happened) THEN
             ! searching from the proc with the least initial cost
             next_in_line: DO ii = 1, nprocs
                acceptable_move = &
                  fb_distribution_acceptable_move(dist(pe(ipe)), &
                                                  dist(pe(ipe))%list(lowest_cost_ind), &
                                                  dist(pe(ii)), &
                                                  average_cost)
                IF ((pe(ii) .NE. pe(ipe)) .AND. acceptable_move) THEN
                   CALL fb_distribution_move(dist(pe(ipe)), &
                                             lowest_cost_ind, &
                                             dist(pe(ii)), &
                                             error)
                   move_happened = .TRUE.
                   EXIT next_in_line
                END IF
             END DO next_in_line
          END IF
          ! if the atom cannot be moved, then this means it is too
          ! costly for all other processes to accept. When this
          ! happens we must stop the redistribution process for this
          ! processor---as all other of its atoms will be even more
          ! costly
          IF (.NOT. move_happened) THEN
             EXIT redistribute
          END IF
       END DO redistribute ! while
    END DO ! ipe

!BEG:DEBUG:LT:2015/06/03
    write (my_pe+200,*) "|LT| after redistribution: "
    DO ipe = 1, nprocs
       write (my_pe+200,*) "|LT| dist(",ipe,") = "
       DO ii = 1, dist(pe(ipe))%nelements
          write (my_pe+200,*) dist(pe(ipe))%list(ii)%id, dist(pe(ipe))%list(ii)%cost
       END DO
       write (my_pe+200,*) "|LT|   cost for this dist = ", dist(pe(ipe))%cost
    END DO
!END:DEBUG:LT:2015/06/03

    DEALLOCATE(pe, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(preferred_procs_set, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! generate local atoms from dist
    ALLOCATE(local_atoms_all(natoms), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(local_atoms_starts(nprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(local_atoms_sizes(nprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_distribution_to_local_atoms(dist, &
                                        local_atoms_all, &
                                        local_atoms_starts, &
                                        local_atoms_sizes, &
                                        error)
    ALLOCATE(local_atoms(local_atoms_sizes(my_pe)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    lb = local_atoms_starts(my_pe)
    ub = local_atoms_starts(my_pe) + local_atoms_sizes(my_pe) - 1
    local_atoms(1:local_atoms_sizes(my_pe)) = local_atoms_all(lb:ub)
    CALL fb_env_set(fb_env=fb_env, &
                    local_atoms=local_atoms, &
                    nlocal_atoms=local_atoms_sizes(my_pe), &
                    error=error)

    ! write out info
    CALL fb_distribution_write_info(dist, scf_section, error)

    DEALLOCATE(local_atoms_all, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(local_atoms_starts, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(local_atoms_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ipe = 1, nprocs
       CALL fb_distribution_release(dist(ipe), error)
    END DO
    DEALLOCATE(dist, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE fb_distribution_build

  PURE FUNCTION fb_distribution_acceptable_move(dist_from, &
                                                element, &
                                                dist_to, &
                                                threshold) &
       RESULT(acceptable)
    TYPE(fb_distribution_list), INTENT(IN) :: dist_from, dist_to
    TYPE(fb_distribution_element), INTENT(IN) :: element
    REAL(KIND=dp), INTENT(IN) :: threshold
    LOGICAL :: acceptable
    
    acceptable = (dist_to%cost + element%cost .LT. dist_from%cost) .AND. &
                 (dist_to%cost .LT. threshold)
  END FUNCTION fb_distribution_acceptable_move

  
  SUBROUTINE fb_distribution_write_info(dist_set, scf_section, error)
    TYPE(fb_distribution_list), DIMENSION(:), INTENT(IN) :: dist_set
    TYPE(section_vals_type), POINTER :: scf_section
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_write_info', &
       routineP = moduleN//':'//routineN
    
    LOGICAL :: failure
    INTEGER :: ii, natoms, nprocs, max_natoms, min_natoms, unit_nr
    REAL(KIND=dp) :: ave_natoms, ave_cost, max_cost, min_cost, total_cost
    TYPE(cp_logger_type), POINTER :: logger

    failure = .FALSE.

    nprocs = SIZE(dist_set)
    natoms = 0
    total_cost = 0.0_dp
    DO ii = 1, nprocs
       natoms = natoms + dist_set(ii)%nelements
       total_cost = total_cost + dist_set(ii)%cost
    END DO
    ave_natoms = REAL(natoms, dp) / REAL(nprocs, dp)
    ave_cost = total_cost / REAL(nprocs, dp)
    max_natoms = 0
    max_cost = 0._dp
    DO ii = 1, nprocs
       max_natoms = MAX(max_natoms, dist_set(ii)%nelements)
       max_cost = MAX(max_cost, dist_set(ii)%cost)
    END DO
    min_natoms = natoms
    min_cost = total_cost
    DO ii = 1, nprocs
       min_natoms = MIN(min_natoms, dist_set(ii)%nelements)
       min_cost = MIN(min_cost, dist_set(ii)%cost)
    END DO
    
    logger => cp_error_get_logger(error=error)
    unit_nr = cp_print_key_unit_nr(logger, scf_section, &
                                   "PRINT%FILTER_MATRIX", &
                                   extension="", &
                                   error=error)

    IF (unit_nr > 0) THEN
       WRITE (UNIT=unit_nr, FMT="(/,A,I6,A)") &
            "FILTER_MAT_DIAG| Load Distribution across ", nprocs, " processors:"
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,A,T55,A,T70,A,T85,A)") &
            "FILTER_MAT_DIAG| ", "Total", "Average", "Max", "Min"
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,I12,T55,F12.1,T70,I12,T85,I10)") &
            "FILTER_MAT_DIAG|   Atomic Matrices", &
            natoms, ave_natoms, max_natoms, min_natoms
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,D12.7,T55,D12.7,T70,D12.7,T85,D12.7)") &
            "FILTER_MAT_DIAG|   Cost*", &
            total_cost, ave_cost, max_cost, min_cost
       WRITE (UNIT=unit_nr, FMT="(A)") &
            "FILTER_MAT_DIAG| (* cost is calculated as sum of cube of atomic matrix sizes)"
    END IF
    CALL cp_print_key_finished_output(unit_nr, logger, scf_section, &
                                      "PRINT%FILTER_MATRIX", error=error)
  END SUBROUTINE fb_distribution_write_info


  SUBROUTINE fb_build_preferred_procs(dbcsr_mat, &
                                      natoms, &
                                      preferred_procs_set, &
                                      common_set_ids, &
                                      n_common_sets, &
                                      error)
    TYPE(cp_dbcsr_type), POINTER :: dbcsr_mat
    INTEGER, INTENT(IN) :: natoms
    TYPE(fb_preferred_procs_list), DIMENSION(:), INTENT(INOUT) :: preferred_procs_set
    INTEGER, DIMENSION(:), INTENT(OUT) :: common_set_ids
    INTEGER, INTENT(OUT) :: n_common_sets
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_build_preferred_procs', &
      routineP = moduleN//':'//routineN

    INTEGER :: icol, pe, prow, nprows, pcol, stat
    LOGICAL :: check_ok, failure
    TYPE(dbcsr_distribution_obj) :: dbcsr_dist
    TYPE(dbcsr_mp_obj) :: dbcsr_mp
    INTEGER, DIMENSION(:), POINTER :: col_dist

    failure = .FALSE.
    check_ok = natoms .LE. cp_dbcsr_nblkcols_total(dbcsr_mat)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    check_ok = SIZE(preferred_procs_set) .GE. natoms
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    check_ok = SIZE(common_set_ids) .GE. natoms
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    CALL cp_dbcsr_get_info(dbcsr_mat, distribution=dbcsr_dist, proc_col_dist=col_dist)
    dbcsr_mp = dbcsr_distribution_mp(dbcsr_dist)
    nprows = dbcsr_mp_nprows(dbcsr_mp)

    DO icol = 1, natoms
       IF (ASSOCIATED(preferred_procs_set(icol)%list)) THEN
          DEALLOCATE(preferred_procs_set(icol)%list, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
       ALLOCATE(preferred_procs_set(icol)%list(nprows), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       pcol = col_dist(icol)
       ! dbcsr prow and pcol counts from 0
       DO prow = 0, nprows-1
          ! here, we count processes from 1, so +1 from mpirank
          preferred_procs_set(icol)%list(prow+1) = &
               dbcsr_mp_get_process(dbcsr_mp, prow, pcol) + 1
       END DO
       preferred_procs_set(icol)%nprocs = nprows
    END DO

    common_set_ids(:) = 0
    common_set_ids(1:natoms) = col_dist(1:natoms) + 1
    n_common_sets = dbcsr_mp_npcols(dbcsr_mp)

  END SUBROUTINE fb_build_preferred_procs


  SUBROUTINE fb_distribution_to_local_atoms(dist_set, &
                                            local_atoms, &
                                            local_atoms_starts, &
                                            local_atoms_sizes, &
                                            error)
    TYPE(fb_distribution_list), DIMENSION(:), INTENT(IN) :: dist_set
    INTEGER, DIMENSION(:), INTENT(OUT) :: local_atoms
    INTEGER, DIMENSION(:), INTENT(OUT) :: local_atoms_starts
    INTEGER, DIMENSION(:), INTENT(OUT) :: local_atoms_sizes
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_to_local_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER :: ipe, iatom, n_procs, pos
    LOGICAL :: check_ok, failure

    failure = .FALSE.

    n_procs = SIZE(dist_set)

    check_ok = SIZE(local_atoms_starts) .GE. n_procs
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    check_ok = SIZE(local_atoms_sizes) .GE. n_procs
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    local_atoms(:) = 0
    local_atoms_starts(:) = 0
    local_atoms_sizes(:) = 0

    pos = 1
    DO ipe = 1, n_procs
       local_atoms_starts(ipe) = pos
       DO iatom = 1, dist_set(ipe)%nelements
          local_atoms(pos) = dist_set(ipe)%list(iatom)%id
          pos = pos + 1
          local_atoms_sizes(ipe) = local_atoms_sizes(ipe) + 1
       END DO
    END DO
  END SUBROUTINE fb_distribution_to_local_atoms


  SUBROUTINE fb_distribution_init(dist, nmax, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    INTEGER, OPTIONAL, INTENT(IN) :: nmax
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_init', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: ii, my_nmax, stat

    failure = .FALSE.
    my_nmax = 0
    IF (PRESENT(nmax)) my_nmax = nmax
    IF (ASSOCIATED(dist%list)) THEN
       DEALLOCATE(dist%list, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    NULLIFY(dist%list)
    IF (my_nmax .GT. 0) THEN
       ALLOCATE(dist%list(my_nmax), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       DO ii = 1, SIZE(dist%list)
          dist%list(ii)%id = 0
          dist%list(ii)%cost = 0.0_dp
       END DO
    END IF
    dist%nelements = 0
    dist%cost = 0.0_dp
  END SUBROUTINE fb_distribution_init


  SUBROUTINE fb_distribution_resize(dist, nmax, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    INTEGER, INTENT(IN) :: nmax
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_resize', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: my_nmax, ii, stat
    TYPE(fb_distribution_element), DIMENSION(:), POINTER :: new_list

    failure = .FALSE.
    IF (.NOT. ASSOCIATED(dist%list)) THEN
       my_nmax = MAX(nmax, 1)
       ALLOCATE(dist%list(my_nmax), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ELSE
       my_nmax = MAX(nmax, dist%nelements)
       ALLOCATE(new_list(my_nmax), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       DO ii = 1, SIZE(new_list)
          new_list(ii)%id = 0
          new_list(ii)%cost = 0.0_dp
       END DO
       DO ii = 1, dist%nelements
          new_list(ii) = dist%list(ii)
       END DO
       DEALLOCATE(dist%list, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       dist%list => new_list
    END IF
  END SUBROUTINE fb_distribution_resize


  SUBROUTINE fb_distribution_add(dist, element, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    TYPE(fb_distribution_element), INTENT(IN) :: element
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_add', &
      routineP = moduleN//':'//routineN

    INTEGER ii, jj, new_nelements, pos
    LOGICAL :: failure

    failure = .FALSE.
    new_nelements = dist%nelements + 1

    ! resize list if necessary
    IF (.NOT. ASSOCIATED(dist%list)) THEN
       CALL fb_distribution_resize(dist, new_nelements, error=error)
    ELSE IF (new_nelements*ENLARGE_RATIO .GT. SIZE(dist%list)) THEN
       CALL fb_distribution_resize(dist, SIZE(dist%list)*EXPAND_FACTOR, error)
    END IF
    ! assuming the list of elements is always sorted with respect to cost
    ! slot the new element into the appropriate spot
    IF (new_nelements == 1) THEN
       dist%list(1) = element
    ELSE
       pos = fb_distribution_find_slot(dist, element)
       DO ii = dist%nelements, pos, -1
          dist%list(ii + 1) = dist%list(ii)
       END DO
       dist%list(pos) = element
    END IF
    dist%nelements = new_nelements
    dist%cost = dist%cost + element%cost
  END SUBROUTINE fb_distribution_add


  PURE FUNCTION fb_distribution_find_slot(dist, element) RESULT(pos)
    TYPE(fb_distribution_list), INTENT(IN) :: dist
    TYPE(fb_distribution_element), INTENT(IN) :: element
    INTEGER :: pos

    INTEGER :: lower, upper, middle, N

    N = dist%nelements
    IF (element%cost .LT. dist%list(1)%cost) THEN
       pos = 1
       RETURN
    END IF
    IF (element%cost .GE. dist%list(N)%cost) THEN
       pos = N + 1
       RETURN
    END IF
    lower = 1
    upper = N
    DO WHILE ((upper - lower) .GT. 1)
       middle = (lower + upper) / 2
       IF (element%cost .LT. dist%list(middle)%cost) THEN
          upper = middle
       ELSE
          lower = middle
       END IF
    END DO
    pos = upper
  END FUNCTION fb_distribution_find_slot


  SUBROUTINE fb_distribution_remove_ind(dist, pos, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    INTEGER, INTENT(IN) :: pos
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_remove_ind', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure, check_ok
    INTEGER :: ii

    failure = .FALSE.
    check_ok = pos .GT. 0
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (pos .LE. dist%nelements) THEN
       dist%cost = dist%cost - dist%list(pos)%cost
       DO ii = pos, dist%nelements-1
          dist%list(ii) = dist%list(ii+1)
       END DO
       dist%list(dist%nelements)%id = 0
       dist%list(dist%nelements)%cost = 0.0_dp
       dist%nelements = dist%nelements - 1
       ! auto resize if required
       IF (dist%nelements*REDUCE_RATIO .LT. SIZE(dist%list)) THEN
          CALL fb_distribution_resize(dist, dist%nelements/SHRINK_FACTOR, error)
       END IF
    END IF
  END SUBROUTINE fb_distribution_remove_ind


  SUBROUTINE fb_distribution_remove_el(dist, element, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    TYPE(fb_distribution_element), INTENT(IN) :: element
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_remove_el', &
      routineP = moduleN//':'//routineN

    INTEGER ii, pos

    pos = dist%nelements + 1
    DO ii = 1, dist%nelements
       IF (element%id == dist%list(ii)%id) THEN
          pos = ii
          EXIT
       END IF
    END DO
    CALL fb_distribution_remove_ind(dist, pos, error)
  END SUBROUTINE fb_distribution_remove_el


  SUBROUTINE fb_distribution_move_ind(dist_from, pos, dist_to, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist_from, dist_to
    INTEGER, INTENT(IN) :: pos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_move_ind', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure
    TYPE(fb_distribution_element) :: element

    failure = .FALSE.
    check_ok = ASSOCIATED(dist_from%list)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    check_ok = pos .LE. dist_from%nelements
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    element = dist_from%list(pos)
    CALL fb_distribution_add(dist_to, element, error)
    CALL fb_distribution_remove(dist_from, pos, error)
  END SUBROUTINE fb_distribution_move_ind


  SUBROUTINE fb_distribution_move_el(dist_from, element, dist_to, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist_from, dist_to
    TYPE(fb_distribution_element), INTENT(IN) :: element
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_move_el', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(dist_from%list)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    CALL fb_distribution_add(dist_to, element, error)
    CALL fb_distribution_remove(dist_from, element, error)
  END SUBROUTINE fb_distribution_move_el


  SUBROUTINE fb_distribution_release(dist, error)
    TYPE(fb_distribution_list), INTENT(INOUT) :: dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_move_el', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (ASSOCIATED(dist%list)) THEN
       DEALLOCATE(dist%list, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
  END SUBROUTINE fb_distribution_release

END Module qs_fb_distribution_methods
