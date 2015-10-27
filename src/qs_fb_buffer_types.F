!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_buffer_types

  USE kinds, ONLY: dp

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_buffer_i_obj, &
            fb_buffer_d_obj

! public methods
  PUBLIC :: fb_buffer_add, &
            fb_buffer_associate, &
            fb_buffer_create, &
            fb_buffer_calc_disps, &
            fb_buffer_calc_sizes, &
            fb_buffer_get, &
            fb_buffer_has_data, &
            fb_buffer_release, &
            fb_buffer_retain, &
            fb_buffer_nullify

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_buffer_types'
  INTEGER, PRIVATE, SAVE :: last_fb_buffer_i_id = 0
  INTEGER, PRIVATE, SAVE :: last_fb_buffer_d_id = 0

! **********************************************************************
!> \brief data for the fb_buffer object (integer)
!> \param n : number of data slices in the buffer
!> \param disps : displacement in data array of each slice, it contains
!>                one more element at the end recording the total
!>                size of the current data, which is the same as the
!>                displacement for the new data to be added
!> \param data : where all of the slices are stored
!> \param id_nr : unique id of this object
!> \param ref_count : reference counter of this object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_buffer_i_data
     INTEGER :: id_nr, ref_count
     INTEGER :: n
     INTEGER, DIMENSION(:), POINTER :: disps
     INTEGER, DIMENSION(:), POINTER :: data
  END type fb_buffer_i_data

! **********************************************************************
!> \brief object/pointer wrapper for fb_buffer object
!> \param obj : pointer to fb_buffer data
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_buffer_i_obj
     TYPE(fb_buffer_i_data), POINTER, PRIVATE :: obj => NULL()
  END type fb_buffer_i_obj

! **********************************************************************
!> \brief data for the fb_buffer object (real, double)
!> \param n : number of data slices in the buffer
!> \param disps : displacement in data array of each slice, it contains
!>                one more element at the end recording the total
!>                size of the current data, which is the same as the
!>                displacement for the new data to be added
!> \param data : where all of the slices are stored
!> \param id_nr : unique id of this object
!> \param ref_count : reference counter of this object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_buffer_d_data
     INTEGER :: id_nr, ref_count
     INTEGER :: n
     INTEGER, DIMENSION(:), POINTER :: disps
     REAL(KIND=dp), DIMENSION(:), POINTER :: data
  END type fb_buffer_d_data

! **********************************************************************
!> \brief object/pointer wrapper for fb_buffer object
!> \param obj : pointer to fb_buffer data
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_buffer_d_obj
     TYPE(fb_buffer_d_data), POINTER, PRIVATE :: obj => NULL()
  END type fb_buffer_d_obj

! method overload interfaces
  INTERFACE fb_buffer_add
     MODULE PROCEDURE fb_buffer_i_add
     MODULE PROCEDURE fb_buffer_d_add
  END INTERFACE fb_buffer_add

  INTERFACE fb_buffer_associate
     MODULE PROCEDURE fb_buffer_i_associate
     MODULE PROCEDURE fb_buffer_d_associate
  END INTERFACE fb_buffer_associate

  INTERFACE fb_buffer_create
     MODULE PROCEDURE fb_buffer_i_create
     MODULE PROCEDURE fb_buffer_d_create
  END INTERFACE fb_buffer_create

  INTERFACE fb_buffer_calc_disps
     MODULE PROCEDURE fb_buffer_i_calc_disps
     MODULE PROCEDURE fb_buffer_d_calc_disps
  END INTERFACE fb_buffer_calc_disps

  INTERFACE fb_buffer_calc_sizes
     MODULE PROCEDURE fb_buffer_i_calc_sizes
     MODULE PROCEDURE fb_buffer_d_calc_sizes
  END INTERFACE fb_buffer_calc_sizes

  INTERFACE fb_buffer_get
     MODULE PROCEDURE fb_buffer_i_get
     MODULE PROCEDURE fb_buffer_d_get
  END INTERFACE fb_buffer_get

  INTERFACE fb_buffer_has_data
     MODULE PROCEDURE fb_buffer_i_has_data
     MODULE PROCEDURE fb_buffer_d_has_data
  END INTERFACE fb_buffer_has_data

  INTERFACE fb_buffer_release
     MODULE PROCEDURE fb_buffer_i_release
     MODULE PROCEDURE fb_buffer_d_release
  END INTERFACE fb_buffer_release

  INTERFACE fb_buffer_retain
     MODULE PROCEDURE fb_buffer_i_retain
     MODULE PROCEDURE fb_buffer_d_retain
  END INTERFACE fb_buffer_retain

  INTERFACE fb_buffer_nullify
     MODULE PROCEDURE fb_buffer_i_nullify
     MODULE PROCEDURE fb_buffer_d_nullify
  END INTERFACE fb_buffer_nullify

CONTAINS

! INTEGER VERSION

! **********************************************************************
!> \brief retains the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_retain(buffer)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_retain', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(buffer%obj))
    buffer%obj%ref_count = buffer%obj%ref_count + 1
  END SUBROUTINE fb_buffer_i_retain

! **********************************************************************
!> \brief releases the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_release(buffer)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(buffer%obj)) THEN
       CPASSERT(buffer%obj%ref_count>0)
       buffer%obj%ref_count = buffer%obj%ref_count - 1
       IF (buffer%obj%ref_count == 0) THEN
          buffer%obj%ref_count = 1
          IF (ASSOCIATED(buffer%obj%data)) THEN
             DEALLOCATE(buffer%obj%data)
          END IF
          IF (ASSOCIATED(buffer%obj%disps)) THEN
             DEALLOCATE(buffer%obj%disps)
          END IF
          buffer%obj%ref_count = 0
          DEALLOCATE(buffer%obj)
       END IF
    ELSE
       NULLIFY(buffer%obj)
    END IF
  END SUBROUTINE fb_buffer_i_release

! *****************************************************************************
!> \brief nullify the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_nullify(buffer)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(buffer%obj)
  END SUBROUTINE fb_buffer_i_nullify

! *****************************************************************************
!> \brief associate object a to object b
!> \param a : object to associate
!> \param b : object target
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_associate(a, b)
    TYPE(fb_buffer_i_obj), INTENT(OUT) :: a
    TYPE(fb_buffer_i_obj), INTENT(IN) :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
    CALL fb_buffer_retain(a)
  END SUBROUTINE fb_buffer_i_associate

! *****************************************************************************
!> \brief check if an object as associated data
!> \param buffer : fb_buffer object
!> \retval res : .TRUE. if buffer has associated data
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  PURE FUNCTION fb_buffer_i_has_data(buffer) RESULT(res)
    TYPE(fb_buffer_i_obj), INTENT(IN) :: buffer
    LOGICAL :: res
    res = ASSOCIATED(buffer%obj)
  END FUNCTION fb_buffer_i_has_data


! *****************************************************************************
!> \brief creates a fb_buffer object
!> \param buffer : fb_buffer object
!> \param max_size : requested total size of the data array
!> \param nslices : total number of slices for the data
!> \param data : the data to be copied to the buffer
!> \param sizes : the size of the slices in the buffer
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_create(buffer,   &
                                max_size, &
                                nslices,  &
                                data,     &
                                sizes)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN), OPTIONAL :: max_size, nslices
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: data, sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_create', &
      routineP = moduleN//':'//routineN

    INTEGER :: my_max_size, my_ndata, my_nslices
    LOGICAL :: check_ok

    ! check optional input
    IF (PRESENT(data)) THEN
       CPASSERT(PRESENT(sizes))
    END IF

    CPASSERT(.NOT.ASSOCIATED(buffer%obj))
    ALLOCATE(buffer%obj)
    ! work out the size of the data array and number of slices
    my_max_size = 0
    my_nslices = 0
    my_ndata = 0
    NULLIFY(buffer%obj%data, &
            buffer%obj%disps)
    ! work out sizes
    IF (PRESENT(max_size)) my_max_size = max_size
    IF (PRESENT(nslices)) my_nslices  = nslices
    IF (PRESENT(sizes)) THEN
       my_nslices = MIN(my_nslices, SIZE(sizes))
       my_ndata = SUM(sizes(1:my_nslices))
       my_max_size = MAX(my_max_size, my_ndata)
    END IF
    ! allocate the arrays
    ALLOCATE(buffer%obj%data(my_max_size))
    ALLOCATE(buffer%obj%disps(my_nslices))
    buffer%obj%data = 0
    buffer%obj%disps = 0
    ! set n for buffer before calc disps 
    buffer%obj%n = my_nslices
    ! compute disps from sizes if required
    IF (PRESENT(sizes)) THEN
       CALL fb_buffer_calc_disps(buffer, sizes)
    END IF
    ! copy data
    IF (PRESENT(data)) THEN
       check_ok = SIZE(data) .GE. my_max_size .AND. &
                  PRESENT(sizes)
       CPASSERT(check_ok)
       buffer%obj%data(1:my_ndata) = data(1:my_ndata)
    END IF
    ! obj meta data update
    buffer%obj%ref_count = 1
    buffer%obj%id_nr = last_fb_buffer_i_id + 1
    last_fb_buffer_i_id = buffer%obj%id_nr
  END SUBROUTINE fb_buffer_i_create

! *****************************************************************************
!> \brief add some data into the buffer
!> \param buffer : fb_buffer object
!> \param data : data to be copied into the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_add(buffer, data)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer
    INTEGER, DIMENSION(:), INTENT(IN) :: data

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_add', &
      routineP = moduleN//':'//routineN

    INTEGER :: new_data_size, new_n, this_size
    INTEGER, DIMENSION(:), POINTER :: new_sizes, new_disps, new_data

    NULLIFY(new_sizes, new_disps, new_data)

    this_size = SIZE(data)
    new_n = buffer%obj%n + 1
    new_data_size = buffer%obj%disps(new_n) + this_size
    ! resize when needed
    IF (SIZE(buffer%obj%disps) .LT. new_n+1) THEN
       ALLOCATE(new_disps(new_n*2))
       new_disps = 0
       new_disps(1:buffer%obj%n+1) = buffer%obj%disps(1:buffer%obj%n+1)
       DEALLOCATE(buffer%obj%disps)
       buffer%obj%disps => new_disps
    END IF
    IF (SIZE(buffer%obj%data) .LT. new_data_size) THEN
       ALLOCATE(new_data(new_data_size*2))
       new_data = 0
       new_data(1:buffer%obj%disps(new_n)) = &
            buffer%obj%data(1:buffer%obj%disps(new_n))
       DEALLOCATE(buffer%obj%data)
       buffer%obj%data => new_data
    END IF
    ! append to the buffer
    buffer%obj%disps(new_n+1) = new_data_size
    buffer%obj%data(buffer%obj%disps(new_n) + 1 : new_data_size) = &
                    data(1:this_size)
    buffer%obj%n = new_n
  END SUBROUTINE fb_buffer_i_add

! *****************************************************************************
!> \brief compute the displacements of each slice in a data buffer from
!>        a given list of sizes of each slice
!> \param buffer : fb_buffer object
!> \param sizes  : list of sizes of each slice on input
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_calc_disps(buffer, sizes)
    TYPE(fb_buffer_i_obj), INTENT(INOUT) :: buffer
    INTEGER, DIMENSION(:), INTENT(IN)    :: sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_calc_disps', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii

    CPASSERT(SIZE(sizes) .GE. buffer%obj%n)
    buffer%obj%disps(1) = 0
    DO ii = 2, buffer%obj%n + 1
       buffer%obj%disps(ii) = buffer%obj%disps(ii-1) + sizes(ii-1)
    END DO
  END SUBROUTINE fb_buffer_i_calc_disps

! *****************************************************************************
!> \brief compute the sizes of each slice
!> \param buffer : fb_buffer object
!> \param sizes  : list of sizes of each slice on output
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_calc_sizes(buffer, sizes)
    TYPE(fb_buffer_i_obj), INTENT(IN)    :: buffer
    INTEGER, DIMENSION(:), INTENT(OUT)   :: sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_calc_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii

    CPASSERT(SIZE(sizes) .GE. buffer%obj%n)
    DO ii = 1, buffer%obj%n
       sizes(ii) = buffer%obj%disps(ii+1) - buffer%obj%disps(ii)
    END DO
  END SUBROUTINE fb_buffer_i_calc_sizes

! *****************************************************************************
!> \brief get data from the fb_buffer object
!> \param buffer  : fb_buffer object
!> \param i_slice : see data
!> \param n     : outputs number of slices in data array
!> \param data_size : outputs the total size of stored data
!> \param sizes : outputs sizes of the slices in data array
!> \param disps : outputs displacements in the data array for each slice
!> \param data  : if i_slice is present:
!>                   returns pointer to the section of data array corresponding
!>                   to i_slice-th slice
!>                else:
!>                   return pointer to the entire non-empty part of the data array
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_i_get(buffer, i_slice, n, data_size, sizes, disps, data)
    TYPE(fb_buffer_i_obj), INTENT(IN)            :: buffer
    INTEGER, INTENT(IN), OPTIONAL                :: i_slice
    INTEGER, INTENT(OUT), OPTIONAL               :: n, data_size
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: sizes, disps
    INTEGER, DIMENSION(:), POINTER, OPTIONAL     :: data

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_get', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(n)) n = buffer%obj%n
    IF (PRESENT(data_size)) data_size = buffer%obj%disps(buffer%obj%n+1)
    IF (PRESENT(sizes)) THEN
       CALL fb_buffer_calc_sizes(buffer, sizes)
    END IF
    IF (PRESENT(disps)) THEN
       CPASSERT(SIZE(disps) .GE. buffer%obj%n)
       disps(1:buffer%obj%n) = buffer%obj%disps(1:buffer%obj%n)
    END IF
    IF (PRESENT(data)) THEN
       IF (PRESENT(i_slice)) THEN
          CPASSERT(i_slice .LE. buffer%obj%n)
          data => buffer%obj%data(buffer%obj%disps(i_slice) + 1 : &
                                  buffer%obj%disps(i_slice+1))
       ELSE
          data => buffer%obj%data(1 : buffer%obj%disps(buffer%obj%n+1))
       END IF
    END IF
  END SUBROUTINE fb_buffer_i_get


! DOUBLE PRECISION VERSION

! **********************************************************************
!> \brief retains the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_retain(buffer)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_retain', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(buffer%obj))
    buffer%obj%ref_count = buffer%obj%ref_count + 1
  END SUBROUTINE fb_buffer_d_retain

! **********************************************************************
!> \brief releases the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_release(buffer)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(buffer%obj)) THEN
       CPASSERT(buffer%obj%ref_count>0)
       buffer%obj%ref_count = buffer%obj%ref_count - 1
       IF (buffer%obj%ref_count == 0) THEN
          buffer%obj%ref_count = 1
          IF (ASSOCIATED(buffer%obj%data)) THEN
             DEALLOCATE(buffer%obj%data)
          END IF
          IF (ASSOCIATED(buffer%obj%disps)) THEN
             DEALLOCATE(buffer%obj%disps)
          END IF
          buffer%obj%ref_count = 0
          DEALLOCATE(buffer%obj)
       END IF
    ELSE
       NULLIFY(buffer%obj)
    END IF
  END SUBROUTINE fb_buffer_d_release

! *****************************************************************************
!> \brief nullify the given fb_buffer
!> \param buffer : the fb_bffer object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_nullify(buffer)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(buffer%obj)
  END SUBROUTINE fb_buffer_d_nullify

! *****************************************************************************
!> \brief associate object a to object b
!> \param a : object to associate
!> \param b : object target
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_associate(a, b)
    TYPE(fb_buffer_d_obj), INTENT(OUT) :: a
    TYPE(fb_buffer_d_obj), INTENT(IN) :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_i_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
    CALL fb_buffer_retain(a)
  END SUBROUTINE fb_buffer_d_associate

! *****************************************************************************
!> \brief check if an object as associated data
!> \param buffer : fb_buffer object
!> \retval res : .TRUE. if buffer has associated data
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  PURE FUNCTION fb_buffer_d_has_data(buffer) RESULT(res)
    TYPE(fb_buffer_d_obj), INTENT(IN) :: buffer
    LOGICAL :: res
    res = ASSOCIATED(buffer%obj)
  END FUNCTION fb_buffer_d_has_data


! *****************************************************************************
!> \brief creates a fb_buffer object
!> \param buffer : fb_buffer object
!> \param max_size : requested total size of the data array
!> \param nslices : total number of slices for the data
!> \param data : the data to be copied to the buffer
!> \param sizes : the size of the slices in the buffer
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_create(buffer,   &
                                max_size, &
                                nslices,  &
                                data,     &
                                sizes)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN), OPTIONAL :: max_size, nslices
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: sizes
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL :: data

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_create', &
      routineP = moduleN//':'//routineN

    INTEGER :: my_max_size, my_ndata, my_nslices
    LOGICAL :: check_ok

    ! check optional input
    IF (PRESENT(data)) THEN
       CPASSERT(PRESENT(sizes))
    END IF

    CPASSERT(.NOT.ASSOCIATED(buffer%obj))
    ALLOCATE(buffer%obj)
    ! work out the size of the data array and number of slices
    my_max_size = 0
    my_nslices = 0
    my_ndata = 0
    NULLIFY(buffer%obj%data, &
            buffer%obj%disps)
    ! work out sizes
    IF (PRESENT(max_size)) my_max_size = max_size
    IF (PRESENT(nslices)) my_nslices  = nslices
    IF (PRESENT(sizes)) THEN
       my_nslices = MIN(my_nslices, SIZE(sizes))
       my_ndata = SUM(sizes(1:my_nslices))
       my_max_size = MAX(my_max_size, my_ndata)
    END IF
    ! allocate the arrays
    ALLOCATE(buffer%obj%data(my_max_size))
    ALLOCATE(buffer%obj%disps(my_nslices+1))
    buffer%obj%data = 0
    buffer%obj%disps = 0
    ! set n for buffer before calc disps
    buffer%obj%n = my_nslices
    ! compute disps from sizes if required
    IF (PRESENT(sizes)) THEN
       CALL fb_buffer_calc_disps(buffer, sizes)
    END IF
    ! copy data
    IF (PRESENT(data)) THEN
       check_ok = SIZE(data) .GE. my_max_size .AND. &
                  PRESENT(sizes)
       CPASSERT(check_ok)
       buffer%obj%data(1:my_ndata) = data(1:my_ndata)
    END IF
    ! obj meta data update
    buffer%obj%ref_count = 1
    buffer%obj%id_nr = last_fb_buffer_i_id + 1
    last_fb_buffer_i_id = buffer%obj%id_nr
  END SUBROUTINE fb_buffer_d_create

! *****************************************************************************
!> \brief add some data into the buffer
!> \param buffer : fb_buffer object
!> \param data : data to be copied into the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_add(buffer, data)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: data

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_add', &
      routineP = moduleN//':'//routineN

    INTEGER :: new_data_size, new_n, this_size
    INTEGER, DIMENSION(:), POINTER :: new_sizes, new_disps
    REAL(KIND=dp), DIMENSION(:), POINTER :: new_data

    NULLIFY(new_sizes, new_disps, new_data)

    this_size = SIZE(data)
    new_n = buffer%obj%n + 1
    new_data_size = buffer%obj%disps(new_n) + this_size
    ! resize when needed
    IF (SIZE(buffer%obj%disps) .LT. new_n+1) THEN
       ALLOCATE(new_disps(new_n*2))
       new_disps = 0
       new_disps(1:buffer%obj%n+1) = buffer%obj%disps(1:buffer%obj%n+1)
       DEALLOCATE(buffer%obj%disps)
       buffer%obj%disps => new_disps
    END IF
    IF (SIZE(buffer%obj%data) .LT. new_data_size) THEN
       ALLOCATE(new_data(new_data_size*2))
       new_data = 0.0_dp
       new_data(1:buffer%obj%disps(new_n)) = &
            buffer%obj%data(1:buffer%obj%disps(new_n))
       DEALLOCATE(buffer%obj%data)
       buffer%obj%data => new_data
    END IF
    ! append to the buffer
    buffer%obj%disps(new_n+1) = new_data_size
    buffer%obj%data(buffer%obj%disps(new_n) + 1 : new_data_size) = &
                    data(1:this_size)
    buffer%obj%n = new_n
  END SUBROUTINE fb_buffer_d_add

! *****************************************************************************
!> \brief compute the displacements of each slice in a data buffer from
!>        a given list of sizes of each slice
!> \param buffer : fb_buffer object
!> \param sizes  : list of sizes of each slice on input
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_calc_disps(buffer, sizes)
    TYPE(fb_buffer_d_obj), INTENT(INOUT) :: buffer
    INTEGER, DIMENSION(:), INTENT(IN)    :: sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_calc_disps', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii

    CPASSERT(SIZE(sizes) .GE. buffer%obj%n)
    buffer%obj%disps(1) = 0
    DO ii = 2, buffer%obj%n + 1
       buffer%obj%disps(ii) = buffer%obj%disps(ii-1) + sizes(ii-1)
    END DO
  END SUBROUTINE fb_buffer_d_calc_disps

! *****************************************************************************
!> \brief compute the sizes of each slice
!> \param buffer : fb_buffer object
!> \param sizes  : list of sizes of each slice on output
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_calc_sizes(buffer, sizes)
    TYPE(fb_buffer_d_obj), INTENT(IN)    :: buffer
    INTEGER, DIMENSION(:), INTENT(OUT)   :: sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_calc_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii

    CPASSERT(SIZE(sizes) .GE. buffer%obj%n)
    DO ii = 1, buffer%obj%n
       sizes(ii) = buffer%obj%disps(ii+1) - buffer%obj%disps(ii)
    END DO
  END SUBROUTINE fb_buffer_d_calc_sizes

! *****************************************************************************
!> \brief get data from the fb_buffer object
!> \param buffer  : fb_buffer object
!> \param i_slice : see data
!> \param n     : outputs number of slices in data array
!> \param data_size : outputs the total size of stored data
!> \param sizes : outputs sizes of the slices in data array
!> \param disps : outputs displacements in the data array for each slice
!> \param data  : if i_slice is present:
!>                   returns pointer to the section of data array corresponding
!>                   to i_slice-th slice
!>                else:
!>                   return pointer to the entire non-empty part of the data array
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_buffer_d_get(buffer, i_slice, n, data_size, sizes, disps, data)
    TYPE(fb_buffer_d_obj), INTENT(IN)            :: buffer
    INTEGER, INTENT(IN), OPTIONAL                :: i_slice
    INTEGER, INTENT(OUT), OPTIONAL               :: n, data_size
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: sizes, disps
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL  :: data

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_buffer_d_get', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(n)) n = buffer%obj%n
    IF (PRESENT(data_size)) data_size = buffer%obj%disps(buffer%obj%n+1)
    IF (PRESENT(sizes)) THEN
       CALL fb_buffer_calc_sizes(buffer, sizes)
    END IF
    IF (PRESENT(disps)) THEN
       CPASSERT(SIZE(disps) .GE. buffer%obj%n)
       disps(1:buffer%obj%n) = buffer%obj%disps(1:buffer%obj%n)
    END IF
    IF (PRESENT(data)) THEN
       IF (PRESENT(i_slice)) THEN
          CPASSERT(i_slice .LE. buffer%obj%n)
          data => buffer%obj%data(buffer%obj%disps(i_slice) + 1 : &
                                  buffer%obj%disps(i_slice+1))
       ELSE
          data => buffer%obj%data(1 : buffer%obj%disps(buffer%obj%n+1))
       END IF
    END IF
  END SUBROUTINE fb_buffer_d_get

END MODULE qs_fb_buffer_types
