!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_trial_fns_types

  USE kinds, ONLY: dp

#include "./base/base_uses.f90"
  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_trial_fns_obj

! public methods
  PUBLIC :: fb_trial_fns_retain,&
            fb_trial_fns_release,&
            fb_trial_fns_nullify,&
            fb_trial_fns_associate,&
            fb_trial_fns_has_data,&
            fb_trial_fns_create,&
            fb_trial_fns_get,&
            fb_trial_fns_set, &
            fb_trial_fns_build_insitu

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_trial_fns_types'
  INTEGER, PRIVATE, SAVE :: last_fb_trial_fns_id = 0

! *****************************************************************************
!> \brief data containing information on trial functions used by filter
!>        matrix diagonalisation method
!> \param nfunctions : nfunctions(ikind) = number of trial functions for
!>                     atomic kind ikind
!> \param functions  : functions(itrial,ikind) = the index of the
!>                     GTO atomic orbital corresponding to itrial-th trial
!>                     function for kind ikind
!> \param id_nr      : unique id for the object
!> \param ref_count  : reference counter for the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_trial_fns_data
     INTEGER :: id_nr, ref_count
     INTEGER, DIMENSION(:), POINTER :: nfunctions
     INTEGER, DIMENSION(:,:), POINTER :: functions
  END TYPE fb_trial_fns_data


! *****************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to fb_trial_fns objects
!> \param obj : pointer to the fb_trial_fns object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_trial_fns_obj
     TYPE(fb_trial_fns_data), POINTER, PRIVATE :: obj
  END TYPE fb_trial_fns_obj

  REAL(KIND=dp), EXTERNAL :: ddot

CONTAINS


! *****************************************************************************
!> \brief retains given object
!> \brief ...
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_retain(trial_fns)
    ! note INTENT(IN) is okay because the obj pointer contained in the
    ! obj type will not be changed
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_retain', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(trial_fns%obj))
    CPASSERT(trial_fns%obj%ref_count>0)
    trial_fns%obj%ref_count = trial_fns%obj%ref_count + 1
  END SUBROUTINE fb_trial_fns_retain


! *****************************************************************************
!> \brief releases given object
!> \brief ...
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_release(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(trial_fns%obj)) THEN
       CPASSERT(trial_fns%obj%ref_count>0)
       trial_fns%obj%ref_count = trial_fns%obj%ref_count - 1
       IF (trial_fns%obj%ref_count == 0) THEN
          trial_fns%obj%ref_count = 1
          IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
             DEALLOCATE(trial_fns%obj%nfunctions)
          END IF
          IF (ASSOCIATED(trial_fns%obj%functions)) THEN
             DEALLOCATE(trial_fns%obj%functions)
          END IF
          trial_fns%obj%ref_count = 0
          DEALLOCATE(trial_fns%obj)
       END IF
    ELSE
       NULLIFY(trial_fns%obj)
    END IF
  END SUBROUTINE fb_trial_fns_release


! *****************************************************************************
!> \brief nullifies the content of given object
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_nullify(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(trial_fns%obj)
  END SUBROUTINE fb_trial_fns_nullify


! *****************************************************************************
!> \brief associates the content of an object to that of another object
!>        of the same type
!> \param a : the output object
!> \param b : the input object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_associate(a, b)
    TYPE(fb_trial_fns_obj), INTENT(OUT)      :: a
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_trial_fns_associate


! *****************************************************************************
!> \brief check if the object has data associated to it
!> \param trial_fns : the fb_trial_fns object in question
!> \retval res : true if trial_fns%obj is associated, false otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_trial_fns_has_data(trial_fns) RESULT(res)
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(trial_fns%obj)
  END FUNCTION fb_trial_fns_has_data


! *****************************************************************************
!> \brief creates an fb_trial_fns object and initialises it
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_create(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_create', &
      routineP = moduleN//':'//routineN

    CPASSERT(.NOT.ASSOCIATED(trial_fns%obj))
    ALLOCATE(trial_fns%obj)
    NULLIFY(trial_fns%obj%nfunctions)
    NULLIFY(trial_fns%obj%functions)
    trial_fns%obj%ref_count = 1
    trial_fns%obj%id_nr = last_fb_trial_fns_id + 1
    last_fb_trial_fns_id = trial_fns%obj%id_nr
  END SUBROUTINE fb_trial_fns_create


! *****************************************************************************
!> \brief initialises an fb_trial_fns object
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_init(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_init', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(trial_fns%obj))
    ! if halo_atoms are associated, then deallocate and de-associate
    IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
       DEALLOCATE(trial_fns%obj%nfunctions)
    END IF
    IF (ASSOCIATED(trial_fns%obj%functions)) THEN
       DEALLOCATE(trial_fns%obj%functions)
    END IF
  END SUBROUTINE fb_trial_fns_init


! *****************************************************************************
!> \brief get values of the attributes of a fb_trial_fns object
!> \param trial_fns  : the fb_trial_fns object in question
!> \param nfunctions : outputs pointer to trial_fns%obj%nfunctions
!> \param functions  : outputs pointer to trial_fns%obj%functions
!> \param max_nfunctions : outputs MAXVAL(nfunctions)
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_get(trial_fns, &
                              nfunctions, &
                              functions, &
                              max_nfunctions)
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nfunctions
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: functions
    INTEGER, OPTIONAL :: max_nfunctions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_get', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(trial_fns%obj))
    IF (PRESENT(nfunctions)) nfunctions => trial_fns%obj%nfunctions
    IF (PRESENT(functions)) functions => trial_fns%obj%functions
    IF (PRESENT(max_nfunctions)) max_nfunctions = MAXVAL(trial_fns%obj%nfunctions)
  END SUBROUTINE fb_trial_fns_get


! *****************************************************************************
!> \brief sets the attributes of a fb_trial_fns object
!> \param trial_fns  : the fb_trial_fns object in question
!> \param nfunctions : associates trial_fns%obj%nfunctions to this pointer
!> \param functions  : associates trial_fns%obj%nfunctions to this pointer
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_set(trial_fns, &
                              nfunctions, &
                              functions)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nfunctions
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: functions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_set', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(trial_fns%obj))
    IF (PRESENT(nfunctions)) THEN
       IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
          DEALLOCATE(trial_fns%obj%nfunctions)
       END IF
       trial_fns%obj%nfunctions => nfunctions
    END IF
    IF (PRESENT(functions)) THEN
       IF (ASSOCIATED(trial_fns%obj%functions)) THEN
          DEALLOCATE(trial_fns%obj%functions)
       END IF
       trial_fns%obj%functions => functions
    END IF
  END SUBROUTINE fb_trial_fns_set


  !!!!!!!!!!!!!!!!!!!

  SUBROUTINE fb_trial_fns_build_insitu(functions, nfunctions, candidate_set)
    INTEGER, DIMENSION(:), INTENT(OUT) :: functions
    INTEGER, INTENT(IN) :: nfunctions
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: candidate_set

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_build_insitu', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: accepted_cols, overlap
    INTEGER :: n_candidates, n_dim, n_accepted, ii
    LOGICAL :: check_ok

    functions = 0
    n_candidates = SIZE(candidate_set, 2)
    n_dim = SIZE(candidate_set, 1)

    ALLOCATE(overlap(n_candidates,n_candidates))
    ALLOCATE(accepted_cols(n_dim, n_candidates))
    overlap = 0.0_dp
    accepted_cols = 0.0_dp

    check_ok = .FALSE.
    n_accepted = 0
    cd: DO ii = 1, n_candidates
       IF (n_accepted == nfunctions) THEN
          EXIT cd
       END IF
       CALL linear_independence_check(candidate_set(:,ii), &
                                      n_dim, &
                                      accepted_cols, &
                                      n_accepted, &
                                      n_candidates, &
                                      overlap, &
                                      check_ok)
       IF (check_ok) THEN
          n_accepted = n_accepted + 1
          functions(n_accepted) = ii
          accepted_cols(:,n_accepted) = candidate_set(:,ii)
       END IF
    END DO cd
    DEALLOCATE(overlap)
    DEALLOCATE(accepted_cols)
    ! check if we have obtained nfunctions
    IF (n_accepted .NE. nfunctions) THEN
       CALL cp_abort(__LOCATION__, &
                     "Unable to find enough trial functions dynamically that"//&
                     "would produce a linearly independent filtered basis set."//&
                     "Consider increase filter termperature or use a larger"//&
                     "basis set for you calculation.")
    END IF
  END SUBROUTINE fb_trial_fns_build_insitu


  SUBROUTINE linear_independence_check(vec, &
                                       dim, &
                                       existing_vecs, &
                                       n_existing, &
                                       n_max, &
                                       overlap, &
                                       res)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vec
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: existing_vecs
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: overlap
    INTEGER, INTENT(IN) :: dim, n_existing, n_max
    LOGICAL, INTENT(OUT) :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'linear_independence_check', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: my_overlap
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work, eigenvalues
    INTEGER :: ii, nn, stat, lwork, info
    REAL(KIND=dp), PARAMETER :: TOLERANCE = 1e-5_dp

    res = .TRUE.

    ! calculating new elements of the overlap matrix in the UPPER triangular
    ! part only. The overlap matrix is symmetric
    nn = n_existing + 1
    DO ii = 1, n_existing
       overlap(ii,nn) = ddot(dim, vec(:), 1, existing_vecs(:,ii), 1)
    END DO
    overlap(nn,nn) = ddot(dim, vec(:), 1, vec(:), 1)

    IF (n_existing == 0) THEN
       RETURN
    END IF

    ! since the number of exisitng_vecs and the dimiension is expected to be
    ! small, use diagonalisation to find eigenvalues here, and we don't care
    ! about work array size, just large enough may be sufficient
    lwork = 3*nn - 1
    ALLOCATE(work(lwork), STAT=stat)
    ! the target matrix will be destroyed affter dsyev, so must make copy
    ALLOCATE(my_overlap(nn,nn), STAT=stat)
    CALL dlacpy("U", nn, nn, overlap, n_max, my_overlap, nn)
    ALLOCATE(eigenvalues(nn), STAT=stat)

    ! the target matrix will be destroyed affter dsyev, so must make copy
    ALLOCATE(my_overlap(nn,nn), STAT=stat)
    CALL dlacpy("U", nn, nn, overlap, n_max, my_overlap, nn)
    ALLOCATE(eigenvalues(nn), STAT=stat)
    ! calculate eigen values
    CALL dsyev("N", "U", nn, my_overlap, nn, eigenvalues, work, lwork, info)

    ! check if all eigenvalues are positive
    DO ii = 1, nn
       IF (eigenvalues(ii) .LT. TOLERANCE) THEN
          res = .FALSE.
          EXIT
       END IF
    END DO

    DEALLOCATE(work, STAT=stat)
    DEALLOCATE(my_overlap, STAT=stat)
    DEALLOCATE(eigenvalues, STAT=stat)
  END SUBROUTINE linear_independence_check

END MODULE qs_fb_trial_fns_types
