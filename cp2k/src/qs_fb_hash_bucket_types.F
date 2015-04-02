MODULE qs_fb_hash_bucket_types

  USE kinds,                           ONLY: int_8

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_hash_bucket_obj

! public methods
  PUBLIC :: fb_hash_bucket_associate, &
            fb_hash_bucket_create, &
            fb_hash_bucket_get, &
            fb_hash_bucket_has_data, &
            fb_hash_bucket_init, &
            fb_hash_bucket_release, &
            fb_hash_bucket_retain, &
            fb_hash_bucket_set

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_hash_bucket_types'
  INTEGER, PRIVATE, SAVE :: last_fb_hash_bucket_id = 0

! *****************************************************************************
!> \brief data defining a bucket in hash table
!> \param key       : key of a hash table entry
!> \param val       : value of a hash table entry
!> \param id_nr     : unique id for the object
!> \param ref_count : reference counter for the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_hash_bucket_data
     INTEGER :: id_nr, ref_count
     INTEGER(KIND=int_8) :: key
     INTEGER :: val
  END TYPE fb_hash_bucket_data


! *****************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to fb_hash_bucket objects
!> \param obj : pointer to the fb_hash_bucket object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_hash_bucket_obj
     TYPE(fb_hash_bucket_data), POINTER, PRIVATE :: obj
  END TYPE fb_hash_bucket_obj


CONTAINS


! *****************************************************************************
!> \brief Associates one fb_hash_bucket object to another
!> \param a : the fb_hash_bucket object to be associated
!> \param b : the fb_hash_bucket object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_associate(a, b)
    TYPE(fb_hash_bucket_obj), INTENT(OUT)    :: a
    TYPE(fb_hash_bucket_obj), INTENT(IN)     :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_hash_bucket_associate


! *****************************************************************************
!> \brief Creates and initialises an empty fb_hash_bucket object
!> \param hash_bucket : the fb_hash_bucket object, its content must be NULL
!>                      and cannot be UNDEFINED
!> \param error       : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_create(hash_bucket, key, val, error)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket
    INTEGER(KIND=int_8), INTENT(IN), &
      OPTIONAL                               :: key
    INTEGER, INTENT(IN), OPTIONAL            :: val
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = .NOT. ASSOCIATED(hash_bucket%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(hash_bucket%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       hash_bucket%obj%key = 0_int_8
       hash_bucket%obj%val = 0
       IF (PRESENT(key)) THEN
          check_ok = PRESENT(val)
          CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
          hash_bucket%obj%key = key
          hash_bucket%obj%val = val
       END IF
       hash_bucket%obj%ref_count = 1
       hash_bucket%obj%id_nr = last_fb_hash_bucket_id + 1
       last_fb_hash_bucket_id = hash_bucket%obj%id_nr
    END IF
  END SUBROUTINE fb_hash_bucket_create


! *****************************************************************************
!> \brief get values of the attributes of a fb_hash_bucket object
!> \param hash_bucket : the fb_hash_bucket object in question
!> \param key         : outputs hash_bucket%obj%key
!> \param val         : outputs hash_bucket%obj%val
!> \param error       : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_get(hash_bucket, &
                                key, &
                                val, &
                                error)
    TYPE(fb_hash_bucket_obj), INTENT(IN)     :: hash_bucket
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: key
    INTEGER, INTENT(OUT), OPTIONAL           :: val
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(hash_bucket%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(key)) key = hash_bucket%obj%key
       IF (PRESENT(val)) val = hash_bucket%obj%val
    END IF
  END SUBROUTINE fb_hash_bucket_get


! *****************************************************************************
!> \brief check if the object has data associated to it
!> \param hash_bucket : the fb_hash_bucket object in question
!> \retval res        : true if hash_bucket%obj is associated, false otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_hash_bucket_has_data(hash_bucket) RESULT(res)
    TYPE(fb_hash_bucket_obj), INTENT(IN)     :: hash_bucket
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(hash_bucket%obj)
  END FUNCTION fb_hash_bucket_has_data


! *****************************************************************************
!> \brief initialises an fb_hash_bucket object
!> \param hash_bucket : the fb_hash_bucket object in question
!> \param error       : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_init(hash_bucket, error)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(hash_bucket%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       hash_bucket%obj%key = 0_int_8
       hash_bucket%obj%val = 0
    END IF
  END SUBROUTINE fb_hash_bucket_init


! *****************************************************************************
!> \brief nullifies the content of given object
!> \param hash_bucket : the fb_hash_bucket object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_nullify(hash_bucket)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(hash_bucket%obj)
  END SUBROUTINE fb_hash_bucket_nullify


! *****************************************************************************
!> \brief releases given object
!> \param hash_bucket : the fb_hash_bucket object in question
!> \param error       : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_release(hash_bucket, error)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(hash_bucket%obj)) THEN
       CPPreconditionNoFail(hash_bucket%obj%ref_count>0, cp_failure_level, routineP, error)
       hash_bucket%obj%ref_count = hash_bucket%obj%ref_count - 1
       IF (hash_bucket%obj%ref_count == 0) THEN
          DEALLOCATE(hash_bucket%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(hash_bucket%obj)
    END IF
  END SUBROUTINE fb_hash_bucket_release


! *****************************************************************************
!> \brief retains given object
!> \param hash_bucket : the fb_hash_bucket object in question
!> \param error       : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_retain(hash_bucket, error)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(hash_bucket%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(hash_bucket%obj%ref_count>0, cp_failure_level, routineP, error)
       hash_bucket%obj%ref_count = hash_bucket%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_hash_bucket_retain


! *****************************************************************************
!> \brief sets the attributes of a fb_hash_bucket object
!> \param hash_bucket  : the fb_hash_bucket object in question
!> \param key          : sets hash_bucket%obj%key to key
!> \param val          : sets hash_bucket%obj%key to val
!> \param error        : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_bucket_set(hash_bucket, &
                                key, &
                                val, &
                                error)
    TYPE(fb_hash_bucket_obj), INTENT(INOUT)  :: hash_bucket
    INTEGER(KIND=int_8), INTENT(IN), &
      OPTIONAL                               :: key
    INTEGER, INTENT(IN), OPTIONAL            :: val
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_bucket_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(hash_bucket%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(key)) hash_bucket%obj%key = key
       IF (PRESENT(val)) hash_bucket%obj%key = val
    END IF
  END SUBROUTINE fb_hash_bucket_set

END MODULE qs_fb_hash_bucket_types
