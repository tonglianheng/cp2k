MODULE qs_fb_env_methods

  USE kinds, ONLY: dp
  USE qs_fb_env_types, ONLY: fb_env_obj, &
                             fb_env_get, &
                             fb_env_set
  USE qs_environment_types,    ONLY: qs_environment_type, &
                                     qs_env_get
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_obj, &
                                     fb_atomic_halo_list_obj, &
                                     fb_atomic_halo_list_init, &
                                     fb_atomic_halo_list_create, &
                                     fb_atomic_halo_list_has_data
  USE particle_types,        ONLY: particle_type
  USE atomic_kind_types,     ONLY: atomic_kind_type, &
                                   get_atomic_kind
  USE qs_neighbor_list,      ONLY: pair_radius_setup
  USE distribution_1d_types, ONLY: distribution_1d_type

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_methods'

  PUBLIC :: fb_env_build_atomic_halos

CONTAINS

  ! **********************************************************************
  !> \brief Builds an fb_atomic_halo_list object using information
  !>        from fb_env
  !> \param fb_env: the fb_env object
  !> \param error: CP2K data container for error logging
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
  SUBROUTINE fb_env_build_atomic_halos(fb_env, error)
    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_atomic_halos', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(qs_environment_type), POINTER :: qs_env
    REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER :: local_particles
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
    INTEGER :: stat, natoms_local, natoms_global, nkinds, ikind, &
               iatom, nkinds_global, jkind, jatom, iatom_local, &
               n_local_atoms_in_kind, n_global_atoms_in_kind, &
               jatom_in_kind, ii, ihalo, ihalo_atom, nhalo_atoms
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    INTEGER, DIMENSION(:), POINTER :: atom_list
    REAL(KIND=dp) :: rij
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: pair_radius
    LOGICAL, DIMENSION(:), ALLOCATABLE :: present_ab

    failure = .FALSE.

    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)

    NULLIFY(qs_env, rcut, particle_set, atomic_kind_set, &
            local_particles, halos, halo_atoms)
    CALL fb_atomic_halo_list_nullify(atomic_halos)

    ! get relevant data from fb_env
    CALL fb_env_get(fb_env=fb_env, &
                    qs_env=qs_env, &
                    rcut=rcut, &
                    atomic_halos=atomic_halos, &
                    error=error)

    IF (fb_atomic_halo_list_has_data(atomic_halos)) THEN
       CALL fb_atomic_halo_list_init(atomic_halos, error)
    ELSE
       CALL fb_atomic_halo_list_create(atomic_halos, error)
    END IF

    ! get the number of atoms and kinds:
    CALL get_qs_env(qs_env=qs_env, &
                    local_particles=local_particles, &
                    natoms=natoms_global, &
                    atomic_kind_set=atomic_kind_set, &
                    nkind=nkinds_global, &
                    error=error)

    nkinds = SIZE(local_particles%list)
    natoms_local = SUM(local_particles%n_el)

!BEG:DEBUG:LT:2014/08/22
    CALL cp_assert(nkinds .LE. SIZE(r_cut), cp_failure_level, &
                   cp_assertion_failed, routineP, &
                   "array size of r_cut and number of local particle kinds mismatch", &
                   error, failure)
    CALL cp_assert(nkinds == nkinds_global, cp_failure_level, &
                   cp_assertion_failed, routineP, &
                   "nkinds not equal to nkinds_global", &
                   error, failure)
!END:DEBUG:LT:2014/08/22

    ! create the halos, one for each local atom
    ALLOCATE(halos(natoms_local), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ihalo = 1, natoms_local
       CALL fb_atomic_halo_create(halos(ihalo), error)
    END DO
    CALL fb_atomic_halo_list_set(atoomic_halos=atomic_halos, &
                                 nhalos=natoms_local, &
                                 error=error)

    ! construct pair_radius
    ALLOCATE(present_ab(nkinds_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(pair_radius(nkinds_global,nkinds_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    present_ab = .TRUE.
    CALL pair_radius_setup(present_ab, present_ab, rcut, rcut, pair_radius, error)
    DEALLOCATE(present_ab, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! calculate the sizes of the halos
    ihalo = 1
    DO ikind=1, nkinds
       n_local_atoms_in_kind = local_particles%n_el(ikind)
       DO iatom_local = 1, n_local_atoms_in_kind
          iatom = local_particles%list(ikind)%array(iatom_local)
          CALL fb_atomic_halo_set(halos(ihalo), &
                                  atom_id=iatom, &
                                  error=error)
          nhalo_atoms = 0
          ! count nhalo_atoms
          DO jkind = 1, nkinds
             CALL get_atomic_kind(atomic_kind=atomic_kind_set(jkind), &
                                  natom=n_global_atoms_in_kind, &
                                  atom_list=atom_list)
             DO jatom_in_kind = 1, n_global_atoms_in_kind
                jatom = atom_list(jatom_in_kind)
                ! calculate distance between iatom and jatom
                rij = 0.0_dp
                DO ii = 1, 3
                   rij = rij + (particle_set(iatom)%r(ii) - particle_set(jatom)%r(ii))**2
                END DO
                rij = SQRT(rij)
                IF (rij .LE. pair_radius(ikind, jkind)) THEN
                   ! jatom is in iatom's halo
                   nhalo_atoms = nhalo_atoms + 1
                END IF
             END DO ! jatom_in_kind
          END DO ! jkind
          CALL fb_atomic_halo_set(halos(ihalo), &
                                  natoms=nhalo_atoms, &
                                  error=error)
          ihalo = ihalo + 1
       END DO ! iatom_local
    END DO ! ikind

    ! once we have the sizes, construct the halos, need to loop over the atoms again
    ihalo = 1
    DO ikind=1, nkinds
       n_local_atoms_in_kind = local_particles%n_el(ikind)
       DO iatom_local = 1, n_local_atoms_in_kind
          iatom = local_particles%list(ikind)%array(iatom_local)
          CALL fb_atomic_halo_get(halos(ihalo), &
                                  natoms=nhalo_atoms, &
                                  error=error)
          ALLOCATE(halo_atoms(nhalo_atoms), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          halo_atoms = 0
          ihalo_atom = 1
          DO jkind = 1, nkinds
             CALL get_atomic_kind(atomic_kind=atomic_kind_set(jkind), &
                                  natom=n_global_atoms_in_kind, &
                                  atom_list=atom_list)
             DO jatom_in_kind = 1, n_global_atoms_in_kind
                jatom = atom_list(jatom_in_kind)
                ! calculate distance between iatom and jatom
                rij = 0.0_dp
                DO ii = 1, 3
                   rij = rij + (particle_set(iatom)%r(ii) - particle_set(jatom)%r(ii))**2
                END DO
                rij = SQRT(rij)
                IF (rij .LE. pair_radius(ikind, jkind)) THEN
                   ! jatom is in iatom's halo
                   halo_atoms(ihalo_atom) = jatom
                   ihalo_atom = ihalo_atom + 1
                END IF
             END DO ! jatom_in_kind
          END DO ! jkind
          CALL fb_atomic_halo_set(halos(ihalo) &
                                  halo_atoms=halo_atoms, &
                                  error=error)
          ihalo = ihalo + 1
       END DO ! iatom_local
    END DO ! ikind

    ! finalise
    CALL fb_atomic_halo_list_set(atomic_halos=atomic_halos, &
                                 halos=halos, &
                                 error=error)

    ! cleanup
    DEALLOCATE(pair_radius, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  END SUBROUTINE fb_atomic_halo_list_build

END MODULE qs_fb_env_methods
