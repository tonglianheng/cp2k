MODULE qs_fb_env_methods

  USE kinds, ONLY: dp
  USE qs_environment_types, ONLY: qs_environment_type
  USE qs_fb_env_types, ONLY: fb_env_type
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_set_type
  USE qs_fb_atomic_halo_methods, ONLY: fb_atomic_halo_set_retain, &
                                       fb_atomic_halo_set_release

#include "cp_common_uses.h"

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_methods'

  INTEGER, PRIVATE, SAVE :: last_fb_env_id = 0

  PUBLIC :: fb_env_retain, &
            fb_env_release, &
            fb_env_create, &
            fb_env_get, &
            fb_env_set

CONTAINS

  ! **********************************************************************
  !> \brief retains the given fb_env
  !> \param fb_env : the fb_env to retain
  !> \param error  : variable to control error logging, stopping, ...
  !>                 see module cp_error_handling
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_retain(fb_env, error)
    TYPE(fb_env_type), POINTER         :: fb_env
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_env_retain', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(fb_env%ref_count>0, cp_failure_level, routineP, error)
       fb_env%ref_count = fb_env%ref_count + 1
    END IF
  END SUBROUTINE fb_env_retain

  ! **********************************************************************
  !> \brief releases a given fb_env
  !> \param fb_env : the fb_env to release
  !> \param error  : variable to control error logging, stopping, ...
  !>                 see module cp_error_handling
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_release(fb_env, error)
    TYPE(fb_env_type), POINTER         :: fb_env
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_env_release', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat

    failure = .FALSE.
    IF (ASSOCIATED(fb_env)) THEN
       CPPreconditionNoFail(fb_env%ref_count>0, cp_failure_level, routineP, error)
       fb_env%ref_count = fb_env%ref_count - 1
       IF (fb_env%ref_count == 0) THEN
          fb_env%ref_count = 1
          IF (ASSOCIATED(fb_env%qs_env)) THEN
             CALL qs_env_release(fb_env%qs_env, error=error)
             NULLIFY(fb_env%qs_env)
          END IF
          IF (ASSOCIATED(fb_env%rcut)) THEN
             DEALLOCATE(fb_env%rcut, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
             NULLIFY(fb_env%rcut)
          END IF
          IF (ASSOCIATED(fb_env%filter_matrix)) THEN
             DEALLOCATE(fb_env%filter_matrix, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
             NULLIFY(fb_env%filter_matrix)
          END IF
          IF (ASSOCIATED(fb_env%atomic_halos)) THEN
             CALL fb_atomic_halo_set_release(fb_env%atomic_halos, error=error)
             NULLIFY(fb_env%atomic_halos)
          END IF
          fb_env%ref_count = 0
          DEALLOCATE(fb_env, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(fb_env)
  END SUBROUTINE fb_env_release


  ! **********************************************************************
  !> \brief creates an empty fb_env_type object
  !> \param fb_env : must be a NULL fb_env_type pointer on input, and output
  !>                 points to an empty fb_env_type object
  !> \param error  : variable to control error logging, stopping, ...
  !>                 see module cp_error_handling
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_create(fb_env, error)
    TYPE(fb_env_type), POINTER         :: fb_env
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER        :: routineN = 'fb_env_create', &
                                          routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    INTEGER                            :: stat
    
    failure = .FALSE.
    CPPrecondition(.NOT. ASSOCIATED(fb_env), cp_failure_level, routineP, error, failure)    
    IF (.NOT. failure) THEN
       ALLOCATE(fb_env, STAT=stat)
       CPPostconditionNoFail(stat==0, cp_failure_level, routineP, error)
       NULLIFY(fb_env%qs_env)
       NULLIFY(fb_env%rcut)
       NULLIFY(fb_env%filter_matrix)
       NULLIFY(fb_env%atomic_halos)
       fb_env%filter_temperature = 0.0_dp
       fb_env%ref_count = 1
       fb_env%id_nr = last_fb_env_id + 1
       last_fb_env_id = fb_env%id_nr
    END IF
  END SUBROUTINE fb_env_create


  ! **********************************************************************
  !> \brief method to get attributes from a given fb_env object
  !> \param fb_env : the fb_env object in question
  !> \param qs_env : outputs pointer to qs_env of fb_env (optional)
  !> \param rcut   : outputs pointer to rcut attribute of fb_env (optional)
  !> \param filter_matrix      : outputs pointer to filter_matrix
  !>                             attribute of fb_env (optional)
  !> \param filter_temperature : outputs pointer to filter_temperature
  !>                             attribute of fb_env (optional)
  !> \param atomic_halos       : outputs pointer to atomic_halos
  !>                             attribute of fb_env (optional)
  !> \param error  : variable to control error logging, stopping, ...
  !>                 see module cp_error_handling
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_get(fb_env, &
                        qs_env, &
                        rcut, &
                        filter_matrix, &
                        filter_temperature, &
                        atomic_halos, &
                        error)
    TYPE(fb_env_type), POINTER :: fb_env
    TYPE(qs_environment_type), POINTER, OPTIONAL :: qs_env
    REAL(kind=dp), DIMENSION(:), POINTER, OPTIONAL :: rcut
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER, OPTIONAL :: filter_matrix
    REAL(kind_dp), OPTIONAL :: filter_temperature
    TYPE(fb_atomic_halo_set_type), POINTER, OPTIONAL :: atomic_halos    
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_get', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(fb_env%ref_count>0, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(qs_env)) qs_env => fb_env%qs_env
       IF (PRESENT(rcut)) rcut => fb_env%rcut
       IF (PRESENT(filter_matrix)) filter_matrix => fb_env%filter_matrix
       IF (PRESENT(filter_temperature)) filter_temperature = fb_env%filter_temperature
       IF (PRESENT(atomic_halos)) atomic_halos => fb_env%atomic_halos
    END IF
  END SUBROUTINE fb_env_get


  ! **********************************************************************
  !> \brief method to set attributes from a given fb_env object
  !> \param fb_env : the fb_env object in question
  !> \param qs_env : sets qs_env attribute of fb_env (optional)
  !> \param rcut   : sets rcut attribute of fb_env (optional)
  !> \param filter_matrix      : sets filter_matrix attribute of fb_env (optional)
  !> \param filter_temperature : sets filter_temperature attribute of fb_env (optional)
  !> \param atomic_halos       : sets atomic_halos attribute of fb_env (optional)
  !> \param error  : variable to control error logging, stopping, ...
  !>                 see module cp_error_handling
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_set(fb_env, &
                        qs_env, &
                        rcut, &
                        filter_matrix, &
                        filter_temperature, &
                        atomic_halos, &
                        error)
    TYPE(fb_env_type), POINTER                   :: fb_env
    TYPE(qs_environment_type), POINTER, OPTIONAL :: qs_env
    TYPE(kind=dp), POINTER,             OPTIONAL :: rcut
    TYPE(cp_dbcsr_p_type), POINTER,     OPTIONAL :: filter_matrix
    REAL(kind_dp), INTENT(in),          OPTIONAL :: filter_temperature
    TYPE(fb_atomic_halo_set_type), POINTER, OPTIONAL :: atomic_halos
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_set', &
                                   routineP = moduleN//':'//routineN

    LOGICAL                            :: failure
    TYPE(cp_error_type), INTENT(inout) :: error

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(fb_env%ref_count>0, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(qs_env)) THEN
          IF (ASSOCIATED(qs_env)) THEN
             CALL qs_env_retain(qs_env, error=error)
          END IF
          CALL qs_env_release(fb_env%qs_env, error=error)
          fb_env%qs_env => qs_env
       END IF
       IF (PRESENT(rcut)) &
            fb_env%rcut => rcut
       IF (PRESENT(filter_matrix)) &
            fb_env%filter_matrix => filter_matrix
       IF (PRESENT(filter_temperature)) &
            fb_env%filter_temperature = filter_temperature
       IF (PRESENT(atomic_halos)) THEN
          IF (ASSOCIATED(atomic_halos)) THEN
             CALL fb_atomic_halo_set_retain(atomic_halos, error=error)
          END IF
          CALL fb_atomic_halo_set_release(fb_env%atomic_halos, error=error)
          fb_env%atomic_halos => atomic_halos
       END IF
    END IF
  END SUBROUTINE fb_env_set
  
END MODULE qs_fb_env_methods
