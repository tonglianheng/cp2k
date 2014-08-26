MODULE qs_fb_atomic_halo_types

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fb_atomic_halo_types'

  PUBLIC :: fb_atomic_halo_type, &
            fb_atomic_halo_p_type, &
            fb_atomic_halo_set_type, &
            fb_atomic_halo_set_p_type

  !**********************************************************************
  !> \brief derived type containing the list of atoms in an atomic halo,
  !>        used by filtered-basis diagonalisation method
  !> \param id_nr      : unique id for the object
  !> \param ref_count  : reference count on the object
  !> \param atom_id    : global atomic id of the atom this halo belongs to
  !> \param natoms     : number of atoms in the halo
  !> \param halo_atoms : the list of global id of atoms in the halo
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  !**********************************************************************
  TYPE fb_atomic_halo_type
     INTEGER :: id_nr, ref_count
     INTEGER :: atom_id
     INTEGER :: natoms
     INTEGER, DIMENSION(:), POINTER :: halo_atoms
  END type fb_atomic_halo_type


  !**********************************************************************
  !> \brief pointer type of fb_atomic_halo_type
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  !**********************************************************************
  TYPE fb_atomic_halo_p_type
     TYPE(fb_atomic_halo_type), POINTER :: obj
  END type fb_atomic_halo_p_type


  !**********************************************************************
  !> \brief derived type describing an atomic halo set used by
  !>        filtered-basis diagonalisation method
  !> \param id_nr     : unique id for the object
  !> \param ref_count : reference count on the object
  !> \param nhalos  : number of halos in the set
  !> \param halos   : halos(ihalo) gives the ihalo-th fb_atomic_halo object
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  !**********************************************************************
  TYPE fb_atomic_halo_set_type
     INTEGER :: id_nr, ref_count
     INTEGER :: nhalos
     TYPE(fb_atomic_halo_p_type), DIMENSION(:), POINTER :: halos
  END type fb_atomic_halo


  !**********************************************************************
  !> \brief pointer type of fb_atomic_halo_set_type
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  !**********************************************************************
  TYPE fb_atomic_halo_set_p_type
     TYPE(fb_atomic_halo_set_type), POINTER :: obj
  END type fb_atomic_halo_set_p_type


END MODULE qs_fb_atomic_halo_types
