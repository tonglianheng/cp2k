MODULE qs_fb_env_types

  USE kinds, ONLY: dp
  USE qs_environment_types, ONLY: qs_environment_type
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_set_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_types'

  PUBLIC :: fb_env_type, &
            fb_env_p_type

  ! **********************************************************************
  !> \brief wrapper to the simulation parameters used for filtered basis
  !>        method
  !> \param qs_env : quickstep environment
  !> \param rcut   : cutoff for included filtered basis set centred at
  !>                 each atom. These defines the ranges of the atomic
  !>                 halos. rcut(ikind) gives the range for atom of
  !>                 global kind ikind
  !> \param filter_matrix : the filter (basis transformation) matrix
  !> \param atomic_halos  : stores information on the neighbors of each
  !>                        atom ii, which are defined by rcut
  !> \param filter_temperature : parameter controlling the smoothness of
  !>                             the filter function during the construction
  !>                             of the filter matrix
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  TYPE fb_env_type
     INTEGER                            :: id_nr, ref_count
     TYPE(qs_environment_type), POINTER :: qs_env
     REAL(kind=dp), DIMENSION(:), POINTER :: rcut
     REAL(kind=dp), DIMENSION(:,:), POINTER :: filter_matrix
     TYPE(fb_atomic_halo_set_type), POINTER :: atomic_halos
     REAL(kind_dp) :: filter_temperature
  END type fb_env_type

  ! *****************************************************************************
  !> \brief allows for the creation of an array of pointers to fb_env
  !> \param obj : pointer to a filtered basis environment
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! *****************************************************************************
  TYPE fb_env_p_type
     TYPE(fb_env_type), POINTER :: obj
  END type fb_env_p_type

END MODULE qs_fb_env_types
