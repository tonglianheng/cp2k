! *****************************************************************************
!> \brief A simple hash table of integer keys, using hash function:
!>          H(k) = (k*p) mod n + 1
!>        where:
!>          k = key
!>          p = a prime number >= n
!>          n = size of the hash table
!>         And collision resolvation is done by open addressing with linear
!>         probing.
!>
!>         The table consists of an array of (key,val) pairs, and
!>         there are no intermediate buckets. For every new entry (k,v):
!>         We first look up slot H(k), and if it already contains an entry,
!>         then move to the next empty slot using a predefined linear probing
!>         sequence (e.g. iterate from slots H(k) to n, and then 1 to H(k)-1).
!>         When we look up, we use the same probing sequence.
! *****************************************************************************
MODULE qs_fb_hash_table_types

  USE kinds, ONLY: int_8

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_hash_table_obj

! public methods
  PUBLIC :: fb_hash_table_add, &
            fb_hash_table_associate, &
            fb_hash_table_create, &
            fb_hash_table_get, &
            fb_hash_table_has_data, &
            fb_hash_table_init, &
            fb_hash_table_rehash, &
            fb_hash_table_release, &
            fb_hash_table_remove, &
            fb_hash_table_retain

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_hash_table_types'
  INTEGER, PRIVATE, SAVE :: last_fb_hash_table_id = 0

! key value indicating an empty slot
  INTEGER(KIND=int_8), PARAMETER, PRIVATE :: EMPTY_KEY = -1_int_8
! Parameters related to automatic resizing of the hash_table:
! Resize by EXPAND_FACTOR if total no. slots / no. of filled slots < ENLARGE_RATIO
  REAL(KIND=dp), PARAMETER, PRIVATE :: ENLARGE_RATIO = 1.0_dp
  REAL(KIND=dp), PARAMETER, PRIVATE :: REDUCE_RATIO = 2.0_dp
  REAL(KIND=dp), PARAMETER, PRIVATE :: EXPAND_FACTOR = 1.5_dp
  REAL(KIND=dp), PARAMETER, PRIVATE :: SHRINK_FACTOR = 1.5_dp

! *****************************************************************************
!> \brief hash table entry data type
!> \param key       : key of the entry
!> \param val       : value of the entry
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_hash_table_element
     INTEGER(KIND=int_8) :: key
     INTEGER :: val
  END type fb_hash_table_element

! *****************************************************************************
!> \brief data defining a hash table using open addressing for collision
!>        resolvation. Uses simple entry structure to be memory efficient
!>        as well as small overhead
!> \param table     : array of buckets of (key, val) pairs
!> \param val       : value of a hash table entry
!> \param id_nr     : unique id for the object
!> \param ref_count : reference counter for the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_hash_table_data
     INTEGER :: id_nr, ref_count
     TYPE(fb_hash_table_element), DIMENSION(:), pointer :: table
     INTEGER :: nelements
     INTEGER :: nmax
     INTEGER :: prime
  END type fb_hash_table_data

! *****************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to fb_hash_table objects
!> \param obj : pointer to the fb_hash_table object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_hash_table_obj
     TYPE(fb_hash_table_data), POINTER, PRIVATE :: obj
  END type fb_hash_table_obj

CONTAINS

! *****************************************************************************
!> \brief Add element to a hash table, auto resize if necessary
!> \param hash_table : the fb_hash_table object
!> \param key        : key of the element
!> \param val        : value of the element
!> \param error      : CP2K error handler container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  RECURSIVE SUBROUTINE fb_hash_table_add(hash_table, key, val, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT) :: hash_table
    INTEGER(KIND=int_8), INTENT(IN) :: key
    INTEGER, INTENT(IN) :: val
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_add', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure
    INTEGER :: islot

    failure = .FALSE.
    check_ok = fb_hash_table_has_data(hash_table)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    ! check hash table size, if too small rehash in a larger table
    IF (INT(hash_table%obj%nelements*ENLARGE_RATIO) .GE. &
        hash_table%obj%nmax) THEN
       CALL fb_hash_table_rehash(hash_table=hash_table, &
                                 nmax=INT(hash_table%obj%nmax*EXPAND_FACTOR), &
                                 error=error)
    END IF
    ! Using open addressing with linear probing to find the first
    ! available empty slot and set (key, val) pair
    islot = fb_hash_table_linear_probe(hash_table, EMPTY_KEY)
    CPPostcondition(islot>0, cp_failure_level, routineP, error, failure)
    hash_table%obj%table(islot)%key = key
    hash_table%obj%table(islot)%val = val
    hash_table%obj%nelements = hash_table%obj%nelements + 1
  END SUBROUTINE fb_hash_table_add

! *****************************************************************************
!> \brief Associates one fb_hash_table object to another
!> \param a : the fb_hash_table object to be associated
!> \param b : the fb_hash_table object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_associate(a, b)
    TYPE(fb_hash_table_obj), INTENT(OUT)    :: a
    TYPE(fb_hash_table_obj), INTENT(IN)     :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_hash_table_associate

! *****************************************************************************
!> \brief Creates and initialises an empty fb_hash_table object
!> \param hash_table : the fb_hash_table object, its content must be NULL
!>                     and cannot be UNDEFINED
!> \param nmax       : total size of the table, optional. If absent default
!>                     size is 1.
!> \param error      : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_create(hash_table, nmax, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)  :: hash_table
    INTEGER, INTENT(IN), OPTIONAL           :: nmax
    TYPE(cp_error_type), INTENT(INOUT)      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_nmax, power, stat
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = .NOT. fb_hash_table_has_data(hash_table)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(hash_table%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       NULLIFY(hash_table%obj%table)
       my_nmax = 1
       IF (PRESENT(nmax)) my_nmax = nmax
       CALL fb_hash_table_init(hash_table=hash_table, &
                               nmax=my_nmax, &
                               error=error)
       ! book keeping stuff
       hash_table%obj%ref_count = 1
       hash_table%obj%id_nr = last_fb_hash_table_id + 1
       last_fb_hash_table_id = hash_table%obj%id_nr
    END IF
  END SUBROUTINE fb_hash_table_create

! *****************************************************************************
!> \brief Retrieve value from a key from a hash table
!> \param hash_table : the fb_hash_table object
!> \param key        : input key
!> \param val        : output value, equals to 0 if key not found
!> \param found      : .TRUE. if key is found, .FALSE. otherwise
!> \param error      : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_get(hash_table, key, val, found, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)  :: hash_table
    INTEGER(KIND=int_8), INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: val
    LOGICAL, INTENT(OUT) :: found
    TYPE(cp_error_type), INTENT(INOUT)      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_get', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure
    INTEGER :: islot

    failure = .FALSE.
    check_ok = fb_hash_table_has_data(hash_table)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    found = .FALSE.
    val = 0
    islot = fb_hash_table_linear_probe(hash_table, key)
    IF (islot > 0) THEN
       val = hash_table%obj%table(islot)%val
       found = .TRUE.
    END IF
  END SUBROUTINE fb_hash_table_get

! *****************************************************************************
!> \brief check if the object has data associated to it
!> \param hash_table : the fb_hash_table object in question
!> \retval res       : true if hash_table%obj is associated, false otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_hash_table_has_data(hash_table) RESULT(res)
    TYPE(fb_hash_table_obj), INTENT(IN)      :: hash_table
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(hash_table%obj)
  END FUNCTION fb_hash_table_has_data

! *****************************************************************************
!> \brief Initialises a fb_hash_table object
!> \param hash_table : the fb_hash_table object, its content must be NULL
!>                     and cannot be UNDEFINED
!> \param nmax       : new size of the table, optional. If absent use the
!>                     old size
!> \param error      : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_init(hash_table, nmax, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)  :: hash_table
    INTEGER, INTENT(IN), OPTIONAL           :: nmax
    TYPE(cp_error_type), INTENT(INOUT)      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_init', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure
    INTEGER :: ii, my_nmax, stat, power

    failure = .FALSE.
    check_ok = fb_hash_table_has_data(hash_table)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    my_nmax = hash_table%obj%nmax
    IF (PRESENT(nmax)) my_nmax = nmax
    ! table length should always be power of 2. Find the least
    ! power that is greater or equal to my_nmax
    power = 0
    DO WHILE (2**power .GE. my_nmax)
       power = power + 1
    END DO
    my_nmax = 2**power
    IF (ASSOCIATED(hash_table%obj%table)) THEN
       IF (SIZE(hash_table%obj%table) .NE. my_nmax) THEN
          DEALLOCATE(hash_table%obj%table, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ALLOCATE(hash_table%obj%table(my_nmax), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       ALLOCATE(hash_table%obj%table(my_nmax), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    hash_table%obj%nmax = my_nmax
    hash_table%obj%prime = fb_hash_table_matching_prime(my_nmax)
    ! initiate element to be "empty"
    DO ii = 1, hash_table%obj%nmax
       hash_table%obj%table(ii)%key = EMPTY_KEY
       hash_table%obj%table(ii)%val = 0
    END DO
    hash_table%obj%nelements = 0
  END SUBROUTINE fb_hash_table_init

! *****************************************************************************
!> \brief Rehash table. If nmax is present, then also change the table size
!>        to MAX(nmax, number_of_non_empty_elements).
!> \param hash_table      : the fb_hash_table object
!> \param nmax [OPTIONAL] : maximum size of the rehashed table
!> \param error           : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_rehash(hash_table, nmax, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)  :: hash_table
    INTEGER, INTENT(IN), OPTIONAL           :: nmax
    TYPE(cp_error_type), INTENT(INOUT)      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_rehash', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: ii, my_nmax, stat
    TYPE(fb_hash_table_element), DIMENSION(:), ALLOCATABLE :: tmp_table

    failure = .FALSE.
    IF (.NOT. fb_hash_table_has_data(hash_table)) THEN
       CALL fb_hash_table_create(hash_table, nmax, error)
       RETURN
    END IF
    IF (PRESENT(nmax)) THEN
       my_nmax = MAX(nmax, hash_table%obj%nelements)
    ELSE
       my_nmax = hash_table%obj%nmax
    END IF
    ALLOCATE(tmp_table(hash_table%obj%nmax), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    tmp_table(:) = hash_table%obj%table(:)
    CALL fb_hash_table_release(hash_table, error)
    CALL fb_hash_table_create(hash_table=hash_table, &
                              nmax=my_nmax, &
                              error=error)
    DO ii = 1, SIZE(tmp_table)
       IF (tmp_table(ii)%key .NE. EMPTY_KEY) THEN
          CALL fb_hash_table_add(hash_table=hash_table, &
                                 key=tmp_table(ii)%key, &
                                 val=tmp_table(ii)%val, &
                                 error=error)
       END IF
    END DO
    DEALLOCATE(tmp_table, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  END SUBROUTINE fb_hash_table_rehash

! *****************************************************************************
!> \brief releases given object
!> \param hash_table : the fb_hash_table object in question
!> \param error      : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_release(hash_table, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)   :: hash_table
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(hash_table%obj)) THEN
       CPPreconditionNoFail(hash_table%obj%ref_count>0, cp_failure_level, routineP, error)
       hash_table%obj%ref_count = hash_table%obj%ref_count - 1
       IF (hash_table%obj%ref_count == 0) THEN
          hash_table%obj%ref_count = 1
          IF (ASSOCIATED(hash_table%obj%table)) THEN
             DEALLOCATE(hash_table%obj%table, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          hash_table%obj%ref_count = 0
          DEALLOCATE(hash_table%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(hash_table%obj)
    END IF
  END SUBROUTINE fb_hash_table_release

! *****************************************************************************
!> \brief Remove element from a table, automatic resize if necessary
!> \param hash_table : the fb_hash_table object
!> \param key        : key of the element to be removed
!> \param error      : CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_remove(hash_table, key, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT) :: hash_table
    INTEGER(KIND=int_8), INTENT(IN) :: key
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_remove', &
      routineP = moduleN//':'//routineN

    LOGICAL :: check_ok, failure
    INTEGER :: islot

    failure = .FALSE.
    check_ok = fb_hash_table_has_data(hash_table)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    islot = fb_hash_table_linear_probe(hash_table, key)
    CPPostcondition(islot>0, cp_failure_level, routineP, error, failure)
    hash_table%obj%table(islot)%key = EMPTY_KEY
    hash_table%obj%nelements = hash_table%obj%nelements - 1
    ! must rehash after setting a filled slot to empty, otherwise the
    ! table will not work. Automatic resize if required
    ! check hash table size, if too small rehash in a larger table
    IF (INT(hash_table%obj%nelements*REDUCE_RATIO) .LT. &
        hash_table%obj%nmax) THEN
       CALL fb_hash_table_rehash(hash_table=hash_table, &
                                 nmax=INT(hash_table%obj%nmax/SHRINK_FACTOR), &
                                 error=error)
    END IF
  END SUBROUTINE fb_hash_table_remove

! *****************************************************************************
!> \brief retains given object
!> \param hash_table : the fb_hash_table object in question
!> \param error      : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_hash_table_retain(hash_table, error)
    TYPE(fb_hash_table_obj), INTENT(INOUT)   :: hash_table
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(hash_table%obj)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(hash_table%obj%ref_count>0, cp_failure_level, routineP, error)
       hash_table%obj%ref_count = hash_table%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_hash_table_retain

! *****************************************************************************
!> \brief Find a prime number equal or larger than ii
!> \param ii   : input integer
!> \retval res : the prime number
! *****************************************************************************
  PURE FUNCTION fb_hash_table_matching_prime(ii) RESULT(res)
    INTEGER, INTENT(IN) :: ii
    INTEGER             :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_hash_table_matching_prime', &
      routineP = moduleN//':'//routineN

    INTEGER :: jj

    ! Find prime by brute force. To check if res is a prime, go
    ! through all integers jj > 2 and < res and see if res mod jj
    ! gives 0. We start from res = ii and then increment by 1 checking
    ! if res is a prime at every step until it is.

    res = ii
    jj = 0
    DO WHILE (jj < res)
       check_prime: DO jj = 2, res-1
          IF (MOD(res, jj) == 0) THEN
             ! res is not a prime, increment
             res = res + 1
             EXIT check_prime
          END IF
       END DO check_prime
    END DO
  END FUNCTION fb_hash_table_matching_prime

! *****************************************************************************
!> \brief Linear probing algorithm for the hash table
!> \param hash_table : the fb_hash_table object
!> \param key        : key to locate
!> \retval islot     : slot location in the table correspond to key, 0 if key
!>                     not found
! *****************************************************************************
  PURE FUNCTION fb_hash_table_linear_probe(hash_table, key) &
       RESULT(islot)
    TYPE(fb_hash_table_obj), INTENT(IN) :: hash_table
    INTEGER(KIND=int_8), INTENT(IN) :: key
    INTEGER :: islot

    INTEGER :: ii, guess

    ! first guess is mapped by the hash_function
    guess = fb_hash_table_hash_function(hash_table, key)
    ! then search for key and stop at first empty slot from guess to
    ! nmax.  using the same linear probe for adding and retreiving
    ! makes all non-empty keys being put before the first empty slot.
    DO ii = guess, hash_table%obj%nmax
       IF (hash_table%obj%table(ii)%key == key) THEN
          islot = ii
          RETURN
       ELSE IF (hash_table%obj%table(ii)%key == EMPTY_KEY) THEN
          islot = 0
          RETURN
       END IF
    END DO
    ! if unsuccessful, search from 1 to guess
    DO ii = 1, guess-1
       IF (hash_table%obj%table(ii)%key == key) THEN
          islot = ii
          RETURN
       ELSE IF (hash_table%obj%table(ii)%key == EMPTY_KEY) THEN
          islot = 0
          RETURN
       END IF
    END DO
    ! if not found set islot to 0
    islot = 0
  END FUNCTION fb_hash_table_linear_probe

! *****************************************************************************
!> \brief Hash function
!> \param hash_table : the fb_hash_table object
!> \param key        : key to locate
!> \retval islot     : slot location in the table correspond to key, 0 if key
!>                     not found
! *****************************************************************************
  PURE FUNCTION fb_hash_table_hash_function(hash_table, key) RESULT(hash)
    TYPE(fb_hash_table_obj), INTENT(IN) :: hash_table
    INTEGER(KIND=int_8), INTENT(IN) :: key
    INTEGER :: hash

    INTEGER(KIND=int_8) :: nmax_8, hash_8, prime_8

    nmax_8 = INT(hash_table%obj%nmax, int_8)
    prime_8 = INT(hash_table%obj%prime, int_8)
    ! IAND here is equivalent to MOD as nmax is alway a power of 2.
    hash_8 = IAND(key*prime_8, nmax_8) + 1_int_8
    hash = INT(hash_8)
  END FUNCTION fb_hash_table_hash_function

END MODULE qs_fb_hash_table_types
