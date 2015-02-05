! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
! *****************************************************************************
  SUBROUTINE  setup_arnoldi_data_s (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_s), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%x_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_s(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_s

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data_s (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER            :: ar_data

    ar_data=>get_data_s(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%x_vec))DEALLOCATE(ar_data%x_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_s

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector_s(arnoldi_data,ind,matrix,vector,error)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj)                          :: vector
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER           :: ar_data
    INTEGER                                           :: vsize, myind, sspace_size, i
    INTEGER, DIMENSION(:), POINTER           :: selected_ind
    COMPLEX(real_4),DIMENSION(:),ALLOCATABLE       :: ritz_v
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    selected_ind=>get_sel_ind(arnoldi_data)
    ar_data=>get_data_s(arnoldi_data)
    sspace_size=get_subsp_size(arnoldi_data)
    vsize=SIZE(ar_data%f_vec)
    myind=selected_ind(ind)
    ALLOCATE(ritz_v(vsize))
    ritz_v=CMPLX(0.0,0.0,real_4)

    IF(dbcsr_is_initialized (vector))CALL dbcsr_release(vector)
    CALL create_col_vec_from_matrix(vector,matrix,1,error)
    IF(control%local_comp)THEN
       DO i=1,sspace_size
          ritz_v(:)=ritz_v(:)+ar_data%local_history(:,i)*ar_data%revec(i,myind)
       END DO
       data_vec => dbcsr_get_data_p (vector%m%data_area, coersion=0.0_real_4)
       ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
       ! again I hope the user knows what is required
       data_vec(1:vsize) =REAL(ritz_v(1:vsize),KIND=real_4)
    END IF

    DEALLOCATE(ritz_v)

  END SUBROUTINE get_selected_ritz_vector_s
     
   
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_initial_vector_s(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_initial_vector_s', &
      routineP = moduleN//':'//routineN
    
    TYPE(arnoldi_data_s), POINTER           :: ar_data
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    INTEGER                                           :: nrow_local, ncol_local
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)

    CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ar_data=>get_data_s(arnoldi_data)
    data_vec => dbcsr_get_data_p (vector%m%data_area, coersion=0.0_real_4)
    IF(nrow_local*ncol_local>0)ar_data%f_vec(1:nrow_local)=data_vec(1:nrow_local)

  END SUBROUTINE set_initial_vector_s  

    
