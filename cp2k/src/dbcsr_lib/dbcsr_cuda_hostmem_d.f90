!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param n size given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_alloc_d (host_mem, n, error)
    REAL(kind=real_8), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL dbcsr_cuda_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*real_8_size, error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    STOP "dbcsr_cuda_hostmem_alloc_d_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE dbcsr_cuda_hostmem_alloc_d



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param n1,n2 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_hostmem_alloc_d_2D (host_mem, n1, n2, error)
    REAL(kind=real_8), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*real_8_size
    CALL dbcsr_cuda_hostmem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    STOP "dbcsr_cuda_hostmem_alloc_d_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE dbcsr_cuda_hostmem_alloc_d_2D


! *****************************************************************************
!> \brief Allocates 3D fortan-array as cuda host-pinned memory.
!> \param n1,n2,n3 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_alloc_d_3D (host_mem, n1, n2, n3, error)
    REAL(kind=real_8), DIMENSION(:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*real_8_size
    CALL dbcsr_cuda_hostmem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3) /))
#else
    STOP "dbcsr_cuda_hostmem_alloc_d_3D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE dbcsr_cuda_hostmem_alloc_d_3D


! *****************************************************************************
!> \brief Allocates 4D fortan-array as cuda host-pinned memory.
!> \param n1,..,n4 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_hostmem_alloc_d_4D (host_mem, n1, n2, n3, n4, error)
    REAL(kind=real_8), DIMENSION(:,:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3, n4
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*MAX(1,n4)*real_8_size
    CALL dbcsr_cuda_hostmem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3),MAX(1,n4) /))
#else
    STOP "dbcsr_cuda_hostmem_alloc_d_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE dbcsr_cuda_hostmem_alloc_d_4D



! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_dealloc_d (host_mem, error)
    REAL(kind=real_8), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_hostmem_dealloc_d', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1)))
    IF (istat /= 0 ) &
       STOP "dbcsr_cuda_hostmem_dealloc_d: Error deallocating host pinned memory"
#else
    STOP "dbcsr_cuda_hostmem_dealloc_d: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_hostmem_dealloc_d


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_2D (host_mem, error)
    REAL(kind=real_8), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_hostmem_dealloc_d_2D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1)))
    IF (istat /= 0 ) &
       STOP "dbcsr_cuda_hostmem_dealloc_d_2D: Error deallocating host pinned memory"
#else
    STOP "dbcsr_cuda_hostmem_dealloc_d: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_2D


! *****************************************************************************
!> \brief Deallocates a 3D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_3D (host_mem, error)
    REAL(kind=real_8), DIMENSION(:,:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_hostmem_dealloc_d_3D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1)))
    IF (istat /= 0 ) &
       STOP "dbcsr_cuda_hostmem_dealloc_d_3D: Error deallocating host pinned memory"
#else
    STOP "dbcsr_cuda_hostmem_dealloc_d: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_3D


! *****************************************************************************
!> \brief Deallocates a 4D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_4D (host_mem, error)
    REAL(kind=real_8), DIMENSION(:,:,:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_hostmem_dealloc_d_4D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1,1)))
    IF (istat /= 0 ) &
       STOP "dbcsr_cuda_hostmem_dealloc_d_4D: Error deallocating host pinned memory"
#else
    STOP "dbcsr_cuda_hostmem_dealloc_d: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_hostmem_dealloc_d_4D

