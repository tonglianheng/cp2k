MODULE qs_fb_atomic_halo_methods

  USE kinds,                   ONLY: dp, &
                                     default_string_length
  USE qs_environment_types,    ONLY: qs_environment_type, &
                                     get_qs_env
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_type, &
                                     fb_atomic_halo_set_type
  USE qs_fb_env_types,         ONLY: fb_env_type, &
                                     fb_env_get
  USE input_section_types,     ONLY: section_vals_type
  USE cp_output_handling,      ONLY: cp_p_file, &
                                     cp_print_key_finished_output, &
                                     cp_print_key_should_output, &
                                     cp_print_key_unit_nr
  USE particle_types,          ONLY: particle_type
  USE atomic_kind_types,       ONLY: atomic_kind_type, &
                                     get_atomic_kind
  USE qs_neighbor_list,        ONLY: pair_radius_setup
  USE distribution_1d_types,   ONLY: distribution_1d_type
  USE cp_para_types,           ONLY: cp_para_env_type
  USE string_utilities,        ONLY: compress

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fb_atomic_halo_methods'

  INTEGER, PRIVATE, SAVE :: last_fb_atomic_halo_id = 0
  INTEGER, PRIVATE, SAVE :: last_fb_atomic_halo_set_id = 0

  PUBLIC :: fb_atomic_halo_retain, &
            fb_atomic_halo_release, &
            fb_atomic_halo_create, &
            fb_atomic_halo_init, &
            fb_atomic_halo_get, &
            fb_atomic_halo_set_retain, &
            fb_atomic_halo_set_release, &
            fb_atomic_halo_set_create, &
            fb_atomic_halo_set_init, &
            fb_atomic_halo_set_build, &
            fb_atomic_halo_set_write

CONTAINS

  SUBROUTINE fb_atomic_halo_retain(atomic_halo, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halo
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_retain', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(atomic_halo%ref_count>0, cp_failure_level, routineP, error)
       atomic_halo%ref_count = fb_env%ref_count + 1
    END IF
  END SUBROUTINE fb_atomic_halo_retain


  SUBROUTINE fb_atomic_halo_release(atomic_halo, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halo
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_release', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat

    failure = .FALSE.
    IF (ASSOCIATED(atomic_halo)) THEN
       CPPreconditionNoFail(atomic_halo%ref_count>0, cp_failure_level, routineP, error)
       atomic_halo%ref_count = atomic_halo%ref_count - 1
       IF (atomic_halo%ref_count == 0) THEN
          atomic_halo%ref_count = 1
          IF (ASSOCIATED(atomic_halo%halo_atoms)) THEN
             DEALLOCATE(atomic_halo%halo_atoms)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halo%ref_count = 0
          DEALLOCATE(atomic_halo, stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(atomic_halo)
  END SUBROUTINE fb_atomic_halo_release


  SUBROUTINE fb_atomic_halo_create(atomic_halo, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halo
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_create', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat

    failure = .FALSE.
    CPPrecondition(.NOT. ASSOCIATED(atomic_halo), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(atomic_halo, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       atomic_halo%atom_id = 0
       atomic_halo%natoms = 0
       NULLIFY(atomic_halo%halo_atoms)
       atomic_halo%ref_count = 1
       atomic_halo%id_nr = last_fb_atomic_halo_id + 1
       last_fb_atomic_halo_id = atomic_halo%id_nr
    END IF
  END SUBROUTINE fb_atomic_halo_create


  SUBROUTINE fb_atomic_halo_init(atomic_halo, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halo
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_init', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ! if halo_atoms are associated, then deallocate and de-associate
       IF (ASSOCIATED(atomic_halo%halo_atoms)) THEN
          DEALLOCATE(atomic_halo%halo_atoms, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(atomic_halo%halo_atoms)
       END IF
       atomic_halo%atom_id = 0
       atomic_halo%natoms = 0
    END IF
  END SUBROUTINE fb_atomic_halo_init


  SUBROUTINE fb_atomic_halo_get(atomic_halo, atom_id, natoms, halo_atoms, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halo
    TYPE(cp_error_type), INTENT(inout) :: error
    INTEGER, INTENT(out), OPTIONAL :: atom_id
    INTEGER, INTENT(out), OPTIONAL :: natoms
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: halo_atoms

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_get', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(atom_id)) atom_id = atomic_halo%atom_id
       IF (PRESENT(natoms)) natoms = atomic_halo%natoms
       IF (PRESENT(halo_atoms)) halo_atoms => atomic_halo%halo_atoms
    END IF
  END SUBROUTINE fb_atomic_halo_get

  SUBROUTINE fb_atomic_halo_set_retain(atomic_halos, error)
    TYPE(fb_atomic_halo_set), POINTER :: atomic_halos
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_set_retain', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(atomic_halos%ref_count>0, cp_failure_level, routineP, error)
       atomic_halo%ref_count = fb_env%ref_count + 1
    END IF
  END SUBROUTINE fb_atomic_halo_set_retain


  SUBROUTINE fb_atomic_halo_set_release(atomic_halos, error)
    TYPE(fb_atomic_halo_type), POINTER :: atomic_halos
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_set_release', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: ii, stat

    failure = .FALSE.
    IF (ASSOCIATED(atomic_halos)) THEN
       CPPreconditionNoFail(atomic_halos%ref_count>0, cp_failure_level, routineP, error)
       atomic_halos%ref_count = atomic_halos%ref_count - 1
       IF (atomic_halos%ref_count == 0) THEN
          atomic_halos%ref_count = 1
          IF (ASSOCIATED(atomic_halos%halos)) THEN
             DO ii = 1, SIZE(halos)
                CALL fb_atomic_halo_release(atomic_halos%halos(ii)%obj, error=error)
             END DO
             DEALLOCATE(atomic_halos%halos, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halos%ref_count = 0
          DEALLOCATE(atomic_halos, stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(atomic_halos)
  END SUBROUTINE fb_atomic_halo_set_release


  SUBROUTINE fb_atomic_halo_set_create(atomic_halos, error)
    TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_set_create', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat

    failure = .FALSE.
    CPPrecondition(.NOT. ASSOCIATED(atomic_halos), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(atomic_halos, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       atomic_halos%nhalos = 0
       NULLIFY(atomic_halos%halos)
       atomic_halos%ref_count = 1
       atomic_halos%id_nr = last_fb_atomic_halo_set_id + 1
       last_fb_atomic_halo_set_id = atomic_halos%id_nr
    END IF
  END SUBROUTINE fb_atomic_halo_set_create


  SUBROUTINE fb_atomic_halo_set_init(atomic_halos, error)
    TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_set_init', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat, ii

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ! if the arrays are associated, then deallocate and de-associate
       IF (ASSOCIATED(atomic_halos%halos)) THEN
          DO ii = 1, SIZE(atomic_halos%halos)
             CALL fb_atomic_halo_release(atomic_halos%halos(ii)%obj, error)
          END DO
          DEALLOCATE(atomic_halos%halos, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(atomic_halo%halos)
       END IF
       atomic_halos%nhalos = 0
    END IF
  END SUBROUTINE fb_atomic_halo_set_init


  SUBROUTINE fb_atomic_halo_set_get(atomic_halos, nhalos, halos, error)
    TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
    INTEGER, INTENT(out), OPTIONAL :: nhalos
    TYPE(fb_atomic_halo_p_type), DIMENSION(:), POINTER, OPTIONAL :: halos
    TYPE(cp_error_type), INTENT(inout) :: error
    
    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_atomic_halo_set_get', &
                                          routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: ihalo, max_nhalo_atoms
    

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(nhalos)) nhalos = atomic_halos%nhalos
       IF (PRESENT(halos)) halos => atomic_halos%halos
    END IF
  END SUBROUTINE fb_atomic_halo_set_get
  

  SUBROUTINE fb_atomic_halo_set_build(atomic_halos, fb_env, error)
    TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
    TYPE(fb_env_type), POINTER :: fb_env
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_set_build', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(qs_environment_type), POINTER :: qs_env
    REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(atomic_halo_type), POINTER :: halo
    TYPE(distribution_1d_type), POINTER :: local_particles
    INTEGER :: stat, natoms_local, natoms_global, nkinds, ikind, &
               iatom, nkinds_global, jkind, jatom, iatom_local, &
               n_local_atoms_in_kind, n_global_atoms_in_kind, &
               jatom_in_kind, ii, ihalo, ihalo_atom
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: pair_radius
    REAL(KIND=dp) :: rij
    LOGICAL, DIMENSION(:), ALLOCATABLE :: present_ab

    failure = .FALSE.

    NULLIFY(qs_env, rcut, particle_set, atomic_kind_set, local_particles, halo)

    ! get relevant data from fb_env
    CALL fb_env_get(fb_env=fb_env, &
                    qs_env=qs_env, &
                    rcut=rcut, &
                    error=error)

    IF (ASSOCIATED(atomic_halos)) THEN
       CALL fb_atomic_halo_set_init(atomic_halos, error)
    ELSE
       CALL fb_atomic_halo_set_create(atomic_halos, error)
    END IF

    ! get the number of atoms and kinds:
    CALL get_qs_env(qs_env=qs_env, &
                    local_particles=local_particles, &
                    natoms=natoms_global, &
                    atomic_kind_set=atomic_kind_set, &
                    nkind=nkinds_global, &
                    error=error)

    nkinds = SIZE(local_particles%list)
    natoms_local = SUM(local_particles%n_el)
!BEG:DEBUG:LT:2014/08/22
    CALL cp_assert(nkinds .LE. SIZE(r_cut), cp_failure_level, &
                   cp_assertion_failed, routineP, &
                   "array size of r_cut and number of local particle kinds mismatch", &
                   error, failure)
    CALL cp_assert(nkinds == nkinds_global, cp_failure_level, &
                   cp_assertion_failed, routineP, &
                   "nkinds not equal to nkinds_global", &
                   error, failure)
!END:DEBUG:LT:2014/08/22

    ! create the halos, one for each local atom
    ALLOCATE(atomic_halos%halos(natoms_local), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ihalo = 1, natoms_local
       CALL fb_atomic_halo_create(atomic_halos%halos(ihalo)%obj, error)
    END DO
    atomic_halos%nhalos = natoms_local

    ! construct pair_radius
    ALLOCATE(present_ab(nkinds_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(pair_radius(nkinds_global,nkinds_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    present_ab = .TRUE.
    CALL pair_radius_setup(present_ab, present_ab, rcut, rcut, pair_radius, error)
    DEALLOCATE(present_ab, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    
    ! calculate the sizes of the halos
    ihalo = 1
    DO ikind=1, nkinds
       n_local_atoms_in_kind = local_particles%n_el(ikind)
       DO iatom_local = 1, n_local_atoms_in_kind
          iatom = local_particles%list(ikind)%array(iatom_local)
          halo => atomic_halos%halos(ihalo)%obj
          halo%atom_id = iatom
          halo%natoms = 0
          DO jkind = 1, nkinds
             CALL get_atomic_kind(atomic_kind=atomic_kind_set(jkind), &
                                  natom=n_global_atoms_in_kind, &
                                  atom_list=atom_list)
             DO jatom_in_kind = 1, n_global_atoms_in_kind
                jatom = atom_list(jatom_in_kind)
                ! calculate distance between iatom and jatom
                rij = 0.0_dp
                DO ii = 1, 3
                   rij = rij + (particle_set(iatom)%r(ii) - particle_set(jatom)%r(ii))**2
                END DO
                rij = SQRT(rij)
                IF (rij .LE. pair_radius(ikind, jkind)) THEN
                   ! jatom is in iatom's halo
                   halo%natoms = halo%natoms + 1
                END IF
             END DO ! jatom_in_kind
          END DO ! jkind
          ihalo = ihalo + 1
       END DO ! iatom_local
    END DO ! ikind

    ! once we have the sizes, construct the halos, need to loop over the atoms again
    DO ihalo = 1, atomic_halos%nhalos
       halo => atomic_halos%halos(ihalo)%obj
       ALLOCATE(halo%halo_atoms(halo%natoms), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END DO
    ihalo = 1
    DO ikind=1, nkinds
       n_local_atoms_in_kind = local_particles%n_el(ikind)
       DO iatom_local = 1, n_local_atoms_in_kind
          iatom = local_particles%list(ikind)%array(iatom_local)
          halo => atomic_halos%halos(ihalo)%obj
          ihalo_atom = 1
          DO jkind = 1, nkinds
             CALL get_atomic_kind(atomic_kind=atomic_kind_set(jkind), &
                                  natom=n_global_atoms_in_kind, &
                                  atom_list=atom_list)
             DO jatom_in_kind = 1, n_global_atoms_in_kind
                jatom = atom_list(jatom_in_kind)
                ! calculate distance between iatom and jatom
                rij = 0.0_dp
                DO ii = 1, 3
                   rij = rij + (particle_set(iatom)%r(ii) - particle_set(jatom)%r(ii))**2
                END DO
                rij = SQRT(rij)
                IF (rij .LE. pair_radius(ikind, jkind)) THEN
                   ! jatom is in iatom's halo
                   halo%halo_atoms(ihalo_atom) = jatom
                   ihalo_atom = ihalo_atom + 1
                END IF
             END DO ! jatom_in_kind
          END DO ! jkind
          ihalo = ihalo + 1
       END DO ! iatom_local
    END DO ! ikind
    
    ! cleanup
    DEALLOCATE(pair_radius, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  END SUBROUTINE fb_atomic_halo_set_build


  SUBROUTINE fb_atomic_halo_set_write(atomic_halos, para_env, fb_section, error)
    TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(section_vals_type), POINTER :: fb_section
    TYPE(cp_error_type), INTENT(inout) :: error
    

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_write', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure, new_file
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: print_unit, mype
    CHARACTER(LEN=default_string_length) :: string
    INTEGER :: ihalo, jatom

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    
    IF (BTEST(cp_print_key_should_output(logger%iter_info, fb_section, &
                                         "PRINT%ATOMIC_HALOS", error=error), &
              cp_p_file)) THEN 
       print_unit = cp_print_key_unit_nr(logger=logger, &
                                         basis_section=fb_section, &
                                         print_key_path="PRINT%ATOMIC_HALOS", &
                                         extension=".out", &
                                         local=.TRUE., &
                                         log_filename=.FALSE., &
                                         file_position="REWIND", &
                                         file_action="WRITE", &
                                         is_new_file=new_file, &
                                         error=error)
       mype = para_env%mepos
       ! print headline
       string = ""
       WRITE (UNIT=string, FMT="(A,I5,A)") &
             "ATOMIC HALOS IN (PROCESS ", mype, ")"
       CALL compress(string)
       IF (pring_unit > 0) THEN
          WRITE (UNIT=print_unit, FMT="(/,/,T2,A)") TRIM(string)
          WRITE (UNIT=print_unit, FMT="(/,T2,A)") &
                "atom : list of atoms in the atomic halo"
       END IF
       ! print content
       DO ihalo = 1, atomic_halos%nhalos
          WRITE (UNIT=print_unit, FMT="(2X,I6,A)", ADVANCE="no") &
                atomic_halos%halos(ihalo)%obj%atom_id, " : "
          DO jatom = 1, atomic_halos%halos(ihalo)%obj%natoms
             WRITE (UNIT=print_unit, FMT="(I6)", ADVANCE="no") &
                   atomic_halos%halos(ihalo)%obj%halo_atoms(jatom)
          END DO
          WRITE (UNIT=print_unit) ""
       END DO
       ! finish
       CALL cp_print_key_finished_output(print_unit, logger, fb_section, &
                                         "PRINT%ATOMIC_HALOS", error=error)
    END IF
  END SUBROUTINE fb_atomic_halo_write

END MODULE qs_fb_atomic_halo_methods
