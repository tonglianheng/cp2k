MODULE qs_fb_matrix_block_types

  USE kinds, ONLY: dp

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: fb_matrix_block_obj

  ! public methods
  PUBLIC :: fb_matrix_block_retain, &
            fb_matrix_block_release, &
            fb_matrix_block_create, &
            fb_matrix_block_has_data, &
            fb_matrix_block_associate, &
            fb_matrix_block_init, &
            fb_matrix_block_get, &
            fb_matrix_block_set

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_matrix_block_types'
  INTEGER, PRIVATE, SAVE :: last_fb_matrix_block_id = 0

  TYPE fb_matrix_block_data
     INTEGER :: id_nr, ref_count
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: el_dp
  END type fb_matrix_block_data
  
  TYPE fb_matrix_block_obj
     TYPE(fb_matrix_block_data), POINTER, PRIVATE :: obj
  END type fb_matrix_block_obj

CONTAINS
  
  SUBROUTINE fb_matrix_block_retain(matrix_block, error)
    TYPE(fb_matrix_block_obj), INTENT(INOUT) :: matrix_block
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix_block%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(matrix_block%obj%ref_count>0, cp_failure_level, routineP, error)
       matrix_block%obj%ref_count = matrix_block%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_matrix_block_retain


  SUBROUTINE fb_matrix_block_release(matrix_block, error)
    TYPE(fb_matrix_block_obj), INTENT(INOUT) :: matrix_block
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(matrix_block%obj)) THEN
       CPPreconditionNoFail(matrix_block%obj%ref_count>0, cp_failure_level, routineP, error)
       matrix_block%obj%ref_count = matrix_block%obj%ref_count - 1
       IF (matrix_block%obj%ref_count == 0) THEN
          matrix_block%obj%ref_count = 1
          IF (ASSOCIATED(matrix_block%obj%el_dp)) THEN
             DEALLOCATE(matrix_block%obj%el_dp, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          matrix_block%obj%ref_count = 0
          DEALLOCATE(matrix_block%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(matrix_block%obj)
    END IF
  END SUBROUTINE fb_matrix_block_release

  
  SUBROUTINE fb_matrix_block_nullify(matrix_block)
    TYPE(fb_matrix_block_obj), INTENT(INOUT)  :: matrix_block

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(matrix_block%obj)
  END SUBROUTINE fb_matrix_block_nullify


  SUBROUTINE fb_matrix_block_associate(a, b)
    TYPE(fb_matrix_block_obj), INTENT(OUT)    :: a
    TYPE(fb_matrix_block_obj), INTENT(IN)     :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_matrix_block_associate


  FUNCTION fb_matrix_block_has_data(matrix_block) RESULT(res)
    TYPE(fb_matrix_block_obj), INTENT(IN)    :: matrix_block
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(matrix_block%obj)
  END FUNCTION fb_matrix_block_has_data


  SUBROUTINE fb_matrix_block_create(matrix_block, input_data, error)
    TYPE(fb_matrix_block_obj), INTENT(INOUT) :: matrix_block
    REAL(KIND=dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: input_data
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_matrix_block_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nrows, ncols, stat
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = .NOT. ASSOCIATED(matrix_block%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(matrix_block%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       NULLIFY(matrix_block%obj%el_dp)
       IF (PRESENT(input_data)) THEN
          nrows = SIZE(input_data, 1)
          ncols = SIZE(input_data, 2)
          ALLOCATE(matrix_block%obj%el_dp(nrows,ncols), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          matrix_block%obj%el_dp(:,:) = input_data(:,:)
       END IF
       matrix_block%obj%ref_count = 1
       matrix_block%obj%id_nr = last_fb_matrix_block_id + 1
       last_fb_matrix_block_id = matrix_block%obj%id_nr
    END IF
  END SUBROUTINE fb_matrix_block_create

  


END MODULE qs_fb_matrix_block_types
