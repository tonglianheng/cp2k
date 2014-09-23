MODULE qs_fb_com_tasks_types
  
  USE kinds, ONLY: int_8

  IMPLICIT NONE
  
#include "cp_common_uses.h"

  PRIVATE

  ! public types
  PUBLIC :: fb_com_tasks_obj

  ! public methods

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_com_tasks_types'
  INTEGER, PRIVATE, SAVE :: last_fb_com_tasks_id = 0

  ! **********************************************************************
  ! explaination on format of task lists (same for tasks_recv and tasks_send):
  ! tasks_recv has dimension (4, ntasks_recv), and stores information on
  ! the block to be copied or transfered
  ! - tasks_recv(TASK_DEST,itask) = destination MPI rank of itask-th task
  ! - tasks_recv(TASK_SRC,itask) = source MPI rank of itask-th task
  ! - tasks_recv(TASK_PAIR,itask) = compressed pair indices of the block of itask-th task
  ! - tasks_recv(TASK_COST,itask) = the cost of itask-th task
  ! 
  ! number of record slots in each task in the task lists
  INTEGER, PARAMETER, PRIVATE :: TASK_N_RECORDS = 4
  ! the indices for the records (1:TASK_DIM) in a task
  INTEGER, PARAMETER, PRIVATE :: TASK_DEST = 1, &
                                 TASK_SRC  = 2, &
                                 TASK_PAIR = 3, &
                                 TASK_COST = 4
  ! **********************************************************************

  TYPE fb_com_tasks_data
     INTEGER :: id_nr, ref_count
     ! use pure integer arrays to facilitate easier MPI coms
     INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks
     INTEGER :: ntasks
  END type fb_com_task_data

  TYPE fb_com_tasks_obj
     TYPE(fb_com_task_data), POINTER, PRIVATE :: obj
  END type fb_com_tasks_obj

CONTAINS

  SUBROUTINE fb_com_tasks_retain(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(in) :: com_tasks
    TYPE(cp_error_type), INTENT(inout) :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_retain', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), &
                   cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(com_tasks%obj%ref_count>0, &
                            cp_failure_level, routineP, error)
       com_tasks%obj%ref_count = com_tasks%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_com_tasks_retain


  SUBROUTINE fb_com_tasks_release(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(inout) :: com_tasks
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_release', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (ASSOCIATED(com_tasks%obj)) THEN
       CPPreconditionNoFail(com_tasks%obj%ref_count>0, &
                            cp_failure_level, routineP, error)
       com_tasks%obj%ref_count = com_tasks%obj%ref_count - 1
       IF (com_tasks%obj%ref_count == 0) THEN
          com_tasks%obj%ref_count = 1
          IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
             DEALLOCATE(com_tasks%obj%tasks)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          com_tasks%obj%ref_count = 0
          DEALLOCATE(com_tasks%obj, stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(com_tasks%obj)
  END SUBROUTINE fb_com_tasks_release


  SUBROUTINE fb_com_tasks_nullify(com_tasks)
    TYPE(fb_com_tasks_obj), INTENT(inout) :: com_tasks

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_nullify', &
                                   routineP = moduleN//':'//routineN

    NULLIFY(com_tasks%obj)
  END SUBROUTINE fb_com_tasks_list_nullify


  SUBROUTINE fb_com_tasks_associate(a, b)
    TYPE(fb_com_tasks_obj), INTENT(in) :: b
    TYPE(fb_com_tasks_obj), INTENT(out) :: a

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_associate', &
                                   routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_com_tasks_associate


  FUNCTION fb_com_tasks_has_data(com_tasks) RESULT(res)
    TYPE(fb_com_tasks_list_obj), INTENT(in) :: com_tasks
    LOGICAL :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_has_data', &
                                   routineP = moduleN//':'//routineN

    res = ASSOCIATED(com_tasks%obj)
  END FUNCTION fb_com_tasks_list_has_data


  
  SUBROUTINE fb_com_tasks_create(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(inout) :: com_tasks
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_create', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(com_tasks%obj), &
                   cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(com_tasks%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       com_tasks%obj%ntasks = 0
       NULLIFY(com_tasks%obj%tasks)
       com_tasks%obj%ref_count = 1
       com_tasks%obj%id_nr = last_fb_com_tasks_id + 1
       last_fb_com_tasks_id = tasks%obj%id_nr
    END IF
  END SUBROUTINE fb_com_tasks_create


  SUBROUTINE fb_com_tasks_init(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(in) :: com_tasks
    TYPE(cp_error_type), INTENT(inout) :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_init', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat
    
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
          DEALLOCATE(com_tasks%obj%tasks, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(com_tasks%obj%tasks)
       END IF
       com_tasks%obj%ntasks = 0
    END IF
  END SUBROUTINE fb_com_tasks_init


  SUBROUTINE fb_com_tasks_get(com_tasks, &
                              ntasks, &
                              tasks, &
                              error)
    TYPE(fb_com_tasks_obj), INTENT(in) :: com_tasks
    INTEGER, INTENT(out) :: ntasks
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks
    TYPE(cp_error_type), INTENT(inout) :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_get', &
                                   routineP = moduleN//':'//routineN
    
    LOGICAL :: failure
    
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(ntasks)) ntasks = com_tasks%obj%ntasks
       IF (PRESENT(tasks)) tasks => com_tasks%obj%tasks
    END IF
  END SUBROUTINE fb_com_tasks_get


  SUBROUTINE fb_com_tasks_set(com_tasks, &
                              ntasks, &
                              tasks, &
                              error)
    TYPE(fb_com_tasks_obj), INTENT(in) :: com_tasks
    INTEGER, INTENT(out) :: ntasks
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks
    TYPE(cp_error_type), INTENT(inout) :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_set', &
                                   routineP = moduleN//':'//routineN
    
    LOGICAL :: failure
    INTEGER :: stat
    
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), &
                   cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(ntasks)) com_tasks%obj%ntasks = ntasks
       IF (PRESENT(tasks)) THEN
          IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
             DEALLOCATE(com_tasks%obj%tasks, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          com_tasks%obj%tasks => tasks
       END IF
    END IF
  END SUBROUTINE fb_com_tasks_set


  SUBROUTINE fb_com_tasks_transpose_dest_src(tasks_dest_is_me, &
                                             n_tasks_dest_is_me, &
                                             direction, &
                                             tasks_src_is_me, &
                                             n_tasks_src_is_me, &
                                             para_env, &
                                             error)
    TYPE(fb_com_tasks_obj), INTENT(inout) :: tasks_dest_is_me
    INTEGER, INTENT(inout) :: n_tasks_dest_is_me
    CHARACTER, INTENT(in) :: direction
    TYPE(fb_com_tasks_obj), INTENT(inout) :: tasks_src_is_me
    INTEGER, INTENT(input) :: n_tasks_src_is_me
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_com_tasks_transpose_dest_src', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: itask, ipe, rank, rank_pos, task_dim, ii, jj, ind, &
               ntasks_in, ntasks_out
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks_in, taks_out
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_sizes, send_disps, &
                                          recv_sizes, recv_disps, &
                                          send_buf, recv_buf
    
    failure = .FALSE.
    
    NULLIFY(tasks_in, tasks_out)
    
    IF (direction == "<") THEN
       CPPrecondition()
       tasks_in => tasks_dest_is_me
       
    ELSE
    END IF


  END SUBROUTINE fb_com_tasks_transpose_dest_src



END MODULE qs_fb_com_tasks_types
