MODULE qs_fb_trial_fns_types
  
  USE kinds, ONLY: dp
  USE qs_integral_utils, ONLY: basis_set_list_setup
  USE qs_kind_types,     ONLY: qs_kind_type
  USE atomic_kind_types, ONLY: atomic_kind_type, &
                               get_atomic_kind
  USE basis_set_types,   ONLY: gto_basis_set_type, &
                               gto_basis_set_p_type, &
                               get_gto_basis_set
  USE orbital_pointers,  ONLY: co, coset

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: fb_trial_fns_obj

  ! public methods
  PUBLIC :: fb_trial_fns_retain, &
            fb_trial_fns_release, &
            fb_trial_fns_nullify, &
            fb_trial_fns_associate, &
            fb_trial_fns_has_data, &
            fb_trial_fns_create, &
            fb_trial_fns_init, &
            fb_trial_fns_get, &
            fb_trial_fns_set
  
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_trial_fns_types'
  INTEGER, PRIVATE, SAVE :: last_fb_trial_fns_id = 0

  TYPE fb_trial_fns_data
     INTEGER :: id_nr, ref_count
     INTEGER, DIMENSION(:), POINTER :: nfunctions
     INTEGER, DIMENSION(:,:), POINTER :: functions
  END type fb_trial_fns_data
  
  TYPE fb_trial_fns_obj
     TYPE(fb_trial_fns_data), POINTER, PRIVATE :: obj
  END type fb_trial_fns_obj


CONTAINS
  
  
  SUBROUTINE fb_trial_fns_retain(trial_fns, error)
    ! note INTENT(in) is okay because the obj pointer contained in the
    ! obj type will not be changed
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_retain', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(trial_fns%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(trial_fns%obj%ref_count>0, cp_failure_level, routineP, error)
       trial_fns%obj%ref_count = trial_fns%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_trial_fns_retain


  SUBROUTINE fb_trial_fns_release(trial_fns, error)
    TYPE(fb_trial_fns_obj), INTENT(inout) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_release', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (ASSOCIATED(trial_fns%obj)) THEN
       CPPreconditionNoFail(trial_fns%obj%ref_count>0, cp_failure_level, routineP, error)
       trial_fns%obj%ref_count = trial_fns%obj%ref_count - 1
       IF (trial_fns%obj%ref_count == 0) THEN
          trial_fns%obj%ref_count = 1
          IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
             DEALLOCATE(trial_fns%obj%nfunctions)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          IF (ASSOCIATED(trial_fns%obj%functions)) THEN
             DEALLOCATE(trial_fns%obj%functions)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          trial_fns%obj%ref_count = 0
          DEALLOCATE(trial_fns%obj, stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(trial_fns%obj)
  END SUBROUTINE fb_trial_fns_release


  
  SUBROUTINE fb_trial_fns_nullify(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(inout) :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_nullify', &
                                   routineP = moduleN//':'//routineN

    NULLIFY(trial_fns%obj)
  END SUBROUTINE fb_trial_fns_nullify



  SUBROUTINE fb_trial_fns_associate(a, b)
    TYPE(fb_trial_fns_obj), INTENT(in) :: b
    TYPE(fb_trial_fns_obj), INTENT(out) :: a

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_associate', &
                                   routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_trial_fns_associate



  FUNCTION fb_trial_fns_has_data(trial_fns) RESULT(res)
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    LOGICAL :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_has_data', &
                                   routineP = moduleN//':'//routineN

    res = ASSOCIATED(trial_fns%obj)
  END FUNCTION fb_trial_fns_has_data



  SUBROUTINE fb_trial_fns_create(trial_fns, error)
    TYPE(fb_trial_fns_obj), INTENT(inout) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_create', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(trial_fns%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(trial_fns%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       NULLIFY(trial_fns%obj%nfunctions)
       NULLIFY(trial_fns%obj%functions)
       trial_fns%obj%ref_count = 1
       trial_fns%obj%id_nr = last_fb_trial_fns_id + 1
       last_fb_trial_fns_id = trial_fns%obj%id_nr
    END IF
  END SUBROUTINE fb_trial_fns_create



  SUBROUTINE fb_trial_fns_init(trial_fns, error)
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_init', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(trial_fns%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ! if halo_atoms are associated, then deallocate and de-associate
       IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
          DEALLOCATE(trial_fns%obj%nfunctions, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(trial_fns%obj%nfunctions)
       END IF
       IF (ASSOCIATED(trial_fns%obj%functions)) THEN
          DEALLOCATE(trial_fns%obj%functions, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(trial_fns%obj%functions)
       END IF
    END IF
  END SUBROUTINE fb_trial_fns_init


  SUBROUTINE fb_trial_fns_get(trial_fns, &
                              nfunctions, &
                              functions, &
                              error)
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: nfunctions
    INTEGER, DIMENSION(:,:), POINTER, OPTIONAL :: functions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_get', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(trial_fns%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(nfunctions)) nfunctions => trial_fns%obj%nfunctions
       IF (PRESENT(functions)) functions => trial_fns%obj%functions
    END IF
  END SUBROUTINE fb_trial_fns_get


  SUBROUTINE fb_trial_fns_set(trial_fns, &
                              nfunctions, &
                              functions, &
                              error)
    TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
    TYPE(cp_error_type), INTENT(inout) :: error
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: nfunctions
    INTEGER, DIMENSION(:,:), POINTER, OPTIONAL :: functions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_set', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(trial_fns%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(nfunctions)) THEN
          IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
             DEALLOCATE(trial_fns%obj%nfunctions, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          trial_fns%obj%nfunctions => nfunctions
       END IF
       IF (PRESENT(functions)) THEN
          IF (ASSOCIATED(trial_fns%obj%functions)) THEN
             DEALLOCATE(trial_fns%obj%functions, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          trial_fns%obj%functions => functions
       END IF
    END IF
  END SUBROUTINE fb_trial_fns_set


  ! SUBROUTINE fb_trial_fns_construct(trial_fns, basis_set_id, qs_kind_set, nelectrons, max_occ, error)
  !   TYPE(fb_trial_fns_obj), INTENT(in) :: trial_fns
  !   INTEGER, INTENT(in) :: basis_set_id
  !   TYPE(qs_kind_type), DIMENSION(:), POINTER :: qs_kind_set
  !   INTEGER, INTENT(in) :: nelectrons
  !   REAL(KIND=dp), INTENT(in) :: max_occ
  !   TYPE(cp_error_type), INTENT(inout) :: error

  !   CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_construct', &
  !                                  routineP = moduleN//':'//routineN

  !   LOGICAL :: failure
  !   TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: basis_set_list
  !   TYPE(gto_basis_set_type), POINTER :: basis_set
  !   INTEGER :: ikind, nkinds, z, icgf, itrial, max_n_trial, counter, stat
  !   INTEGER, DIMENSION(:), POINTER :: nfunctions
  !   INTEGER, DIMENSION(:,:), POINTER :: functions
  !   REAL(KIND=dp), POINTER :: nelec_float

  !   NULLIFY(nfunctions)

  !   nkinds = SIZE(qs_kind_set)
  !   nelec_float = REAL(nelectrons, dp)

  !   ALLOCATE(nfunctions(nkinds), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  !   ALLOCATE(basis_set_list(nkinds), STAT=istat)
  !   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  !   CALL basis_set_list_setup(basis_set_list, basis_set_id, qs_kind_set)

  !   DO ikind = 1, nkinds
  !      CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
  !                           z=z)
  !      ! "gto = gaussian type orbital"
  !      basis_set => basis_set_list(ikind)%gto_basis_set
  !      CALL get_gto_basis_set(gto_basis_set=basis_set, &
  !                             nset=nset, &
  !                             lmax=lmax)
  !      basis_set: DO iset = 1, nset
  !         old_lshell = lmax
  !         DO ishell = 1, nshell(iset)
  !            lshell = gto_basis_set%l(ishell,iset)
  !            counter = 0
  !            ! loop over orbitals within the same l
  !            DO ico = ncoset(lshell-1)+1, ncoset(lshell)
  !               counter = counter + 1
  !               ! only include the first zeta orbitals
  !               IF (counter .LE. nco(lshell)) THEN
  !                  nfunctions(ikind) = nfunctions(ikind) + 1
  !               END IF
  !            END DO
  !            ! we have got enough trial functions when we have more
  !            ! basis functions than the number of electrons (obtained
  !            ! from atomic z), AND that that we have included all the
  !            ! first zeta orbitals of an angular momentum quantum
  !            ! number l
  !            IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax)) .AND. &
  !                (max_occ * nfunctions(ikind) .GE. nelec_float)) THEN
  !               EXIT basis_set
  !            END IF
  !            old_lshell = lshell
  !         END DO
  !      END DO basis_set
  !   END DO ! ikind  

  !   ! now that we have the number of trial functions get the trial
  !   ! functions
  !   max_n_trial = MAXVAL(nfunctions)
  !   ALLOCATE(functions(max_n_trial,nkinds), STAT=stat)
  !   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  !   functions(:,:) = 0
  !   ! redo the loops to get the trial function indices within the basis set
  !   DO ikind = 1, nkinds
  !      CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
  !                           z=z)
  !      ! "gto = gaussian type orbital"
  !      basis_set => basis_set_list(ikind)%gto_basis_set
  !      CALL get_gto_basis_set(gto_basis_set=gto_basis_set, &
  !                             nset=nset, &
  !                             lmax=lmax)
  !      icgf = 0
  !      basis_set: DO iset = 1, nset
  !         old_lshell = lmax
  !         DO ishell = 1, nshell(iset)
  !            lshell = gto_basis_set%l(ishell,iset)
  !            counter = 0
  !            ! loop over orbitals within the same l
  !            DO ico = ncoset(lshell-1)+1, ncoset(lshell)
  !               icgf = icgf + 1
  !               counter = counter + 1
  !               ! only include the first zeta orbitals
  !               IF (counter .LE. nco(lshell)) THEN
  !                  itrial = itrial + 1
  !                  functions(itrial,ikind) = icgf
  !               END IF
  !            END DO
  !            ! we have got enough trial functions when we have more
  !            ! basis functions than the number of electrons (obtained
  !            ! from atomic z), AND that that we have included all the
  !            ! first zeta orbitals of an angular momentum quantum
  !            ! number l
  !            IF (((lshell .GT. old_lshell) .OR. &
  !                 (lshell .EQ. lmax)) .AND. &
  !                (nfunctions(ikind) .GE. z)) THEN
  !               EXIT basis_set
  !            END IF
  !            old_lshell = lshell
  !         END DO
  !      END DO basis_set
  !   END DO ! ikind  

  !   ! set trial_functions
  !   CALL fb_trial_fns_set(trial_fns=trial_fns, &
  !                         nfunctions=nfunctions, &
  !                         functions=functions, &
  !                         error=error)

  !   ! cleanup
  !   DEALLOCATE(basis_set_list, STAT=stat)
  !   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  ! END SUBROUTINE fb_trial_fns_construct


  ! FUNCTION fb_trial_fns_calc_gto_radius(gto_basis_set, l, threshold, epsin) RESULT(radius)

  !   DO ipgf = npgf
  !      radius = MAX(radius, exp_radius(l, alpha(ipgf), threshold, contract_coeff(ipgf), espin))
  !   END DO
    
  ! END FUNCTION fb_trial_fns_calc_gto_radius

END MODULE qs_fb_trial_fns_types
