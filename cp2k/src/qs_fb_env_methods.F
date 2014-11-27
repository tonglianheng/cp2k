MODULE qs_fb_env_methods

  USE kinds, ONLY: dp, &
                   int_8
  USE qs_fb_env_types, ONLY: fb_env_obj, &
                             fb_env_get, &
                             fb_env_set
  USE qs_environment_types,    ONLY: qs_environment_type, &
                                     qs_env_get
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_obj, &
                                     fb_atomic_halo_get, &
                                     fb_atomic_halo_set, &
                                     fb_atomic_halo_sort, &
                                     fb_atomic_halo_nelectrons_estimate_Z, &
                                     fb_atomic_halo_list_obj, &
                                     fb_atomic_halo_list_init, &
                                     fb_atomic_halo_list_create, &
                                     fb_atomic_halo_list_get, &
                                     fb_atomic_halo_list_set, &
                                     fb_atomic_halo_list_has_data
  USE particle_types,        ONLY: particle_type
  USE atomic_kind_types,     ONLY: atomic_kind_type, &
                                   get_atomic_kind
  USE qs_neighbor_list,      ONLY: pair_radius_setup
  USE distribution_1d_types, ONLY: distribution_1d_type
  USE cp_dbcsr_types, ONLY: cp_dbcsr_type
  USE cp_dbcsr_interface, ONLY: cp_dbcsr_get_stored_coordinates, &
                                cp_dbcsr_get_info
  USE dbcsr_util, ONLY: convert_sizes_to_offsets
  USE array_types, ONLY: array_data
  USE cp_para_types, ONLY: cp_para_env_type
  USE task_list_types, ONLY: task_list_type
  USE memory_utilities, ONLY: reallocate
  USE fermi_utils, ONLY: Fermi, &
                         FermiFixed
  USE qs_mo_types, ONLY: mo_set_p_type, &
                         get_mo_set

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
  SUBROUTINE fb_env_build_atomic_halos(fb_env, qs_env, error)
    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_atomic_halos', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER :: local_particles
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
    INTEGER :: stat, natoms_local, natoms_global, nkinds, ikind, &
               iatom, nkinds_global, jkind, jatom, iatom_local, &
               n_local_atoms_in_kind, n_global_atoms_in_kind, &
               jatom_in_kind, ii, ihalo, ihalo_atom, nhalo_atoms, &
               nelectrons_in_halo, owner_id_in_halo
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
                    rcut=rcut, &
                    atomic_halos=atomic_halos, &
                    error=error)

    ! create atomic_halos
    CALL fb_atomic_halo_list_create(atomic_halos, error)

    ! get the number of atoms and kinds:
    CALL get_qs_env(qs_env=qs_env, &
                    local_particles=local_particles, &
                    natoms=natoms_global, &
                    particle_set=particle_set, &
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

    ! build halos
    ihalo = 1
    DO ikind=1, nkinds
       n_local_atoms_in_kind = local_particles%n_el(ikind)
       DO iatom_local = 1, n_local_atoms_in_kind
          iatom = local_particles%list(ikind)%array(iatom_local)
          CALL fb_atomic_halo_set(halos(ihalo), &
                                  owner_atom=iatom, &
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
          ! allocate
          ALLOCATE(halo_atoms(nhalo_atoms), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ! now find and store halo atoms
          ihalo_atom = 1
          DO jkind = 1, nkinds
             CALL get_atomic_kind(atomic_kind=atomic_kind_set(jkind), &
                                  natom=n_global_atoms_in_kind, &
                                  atom_list=atom_list)
             DO jatom_in_kind = 1, n_global_atoms_in_kind
                jatom = atom_list(jatom_in_kind)
                ! record the position of the atom owning the halo in its own halo
                IF (jatom == iatom) THEN
                   owner_id_in_halo = ihalo_atom
                END IF
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
          ! calculate the number of electrons in each halo
          nelectrons_in_halo = fb_atomic_halo_nelectrons_estimate_Z(halos(ihalo), &
                                                                    particle_set) 
          ! set atomic halo
          CALL fb_atomic_halo_set(halos(ihalo) &
                                  owner_id_in_halo=owner_id_in_halo, &
                                  natoms=nhalo_atoms, &
                                  halo_atoms=halo_atoms, &
                                  nelectrons=nelectrons_in_halo, &
                                  error=error)
          ! sort atomic halo
          CALL fb_atomic_halo_sort(halos(ihalo), error)XS

          ihalo = ihalo + 1
       END DO ! iatom_local
    END DO ! ikind

    ! finalise
    CALL fb_atomic_halo_list_set(atomic_halos=atomic_halos, &
                                 halos=halos, &
                                 error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    atomic_halos=atomic_halos, &
                    error=error)

    ! cleanup
    DEALLOCATE(pair_radius, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  END SUBROUTINE fb_env_build_atomic_halos


  SUBROUTINE fb_env_build_filter_matrix(fb_env, qs_env, error)
    TYPE(fb_env_obj), POINTER :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_env_build_filter_matrix', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure

    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(dbcsr_distribution_obj), POINTER :: dbcsr_dist
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mo_list
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: mat_S, mat_H
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: atomic_H, atomic_S
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: filtered_eigenvectors, &
                                                  atomic_filter_matrix
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: filter_function, &
                                                eigenvalues
    REAL(KIND=dp) :: filter_temperature, fermi_level, maxocc
    INTEGER :: ispin, nblkrows_total, nblkcols_total, stat, jhalo, nhalos, &
               this_atom, this_atom_in_halo, ii, jj, mat_dim, nrows, ncols, &
               natoms_in_halo, ihalo_atom, jhalo_atom, iatom_global, jatom_global, &
               natoms_global, ikind, jkind, nspin
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    INTEGER, DIMENSION(:), ALLOCATABLE :: atomic_H_blk_row_start, &
                                          atomic_H_blk_col_start, &
                                          atomic_S_blk_row_start, &
                                          atomic_S_blk_col_start
    INTEGER, DIMENSION(:), POINTER :: cbs
    TYPE(array_i1d_obj) :: col_blk_size, row_blk_size

    failure = .FALSE.

    CPPostcondition(fb_env_has_data(fb_env), &
                    cp_failure_level, routineP, error, failure)

    NULLIFY(para_env, mo_list, particle_set, atomic_kind, mat_S, mat_H, halos, &
            atomic_H, atomic_S, atomic_H_blk_row_start, atomic_H_blk_col_start, &
            atomic_S_blk_row_start, atomic_S_blk_col_start, cbs)
    CALL array_nullify(row_blk_size)
    CALL array_nullify(col_blk_size)
    CALL fb_atomic_halo_list_nullify(atomic_halos)

    CALL fb_env_get(fb_env=fb_env, &
                    filter_temperature=filter_temperature, &
                    error=error)

    CALL get_qs_env(qs_env=qs_env, &
                    para_env=para_env, &
                    mos=mo_list, &
                    particle_set=particle_set, &
                    matrix_s=mat_S, &
                    matrix_h=mat_H, &
                    dbcsr_dist=dbcsr_dist, &
                    error=error)

    CPPrecondition(SIZE(mat_H)==SIZE(mat_S), cp_failure_level, routineP, error, failure)
    nspin = SIZE(mat_H)

    ! get trial function information
    CALL fb_trial_fns_get(trial_fns=trial_fns, &
                          nfns=nfns, &
                          fns=fns, &
                          error=error)
    
    ! allocate the row_blk_sizes and col_blk_sizes arrays for
    ! constructing the filter matrix in DBCSR format
    natoms_global = SIZE(particle_set)
    ALLOCATE(cbs(natoms_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    rbs = 0
    cbs = 0
    ! row_blk_size is same as S or H matrix
    CALL cp_dbcsr_get_info(mat_S(ispin=1)%matrix, &
                           nblkrows_total=nblkrows_total, &
                           nblkcols_total=nblkcols_total, &
                           row_blk_size=row_blk_size)
    ! work out col_blk_size
    DO iatom_global = 1, natoms_global
       atomic_kind => particle_set(iatom_global)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=ikind)
       cbs(iatom_global) = nfns(ikind)
    END DO
    ! do not deallocate cbs if gift=.TRUE. as col_blk_sizes will only point to cbs
    CALL array_new(col_blk_sizes, cbs, gift=.TRUE.)
    ! the filter matrix is non-square and is always non-symmetric
    symmetry_string = dbcsr_type_no_symmetry
    ! create empty filter matrices for each spin
    DO ispin = 1, nspin
       WRITE (spin_string, FMT="(I1)"), ispin
       name = TRIM("FILTER MATRIX SPIN"//spin_string)
       CALL compress(name)
       CALL uppercase(name)
       ALLOCATE(filter_matrix(ispin)%matrix, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       CALL cp_dbcsr_init(filter_matrix(ispin)%matrix, error=error)
       CALL cp_dbcsr_create(matrix=filter_matrix(ispin)%matrix, &
                            name=name, &
                            dist=dbcsr_dist, &
                            matrix_type=symmetry_string, &
                            row_blk_size=row_blk_size, &
                            col_blk_size=col_blk_size, &
                            nblks=0, &
                            nze=0, &
                            error=error)
    END DO
    ! generate atomic halos
    CALL fb_env_build_atomic_halos(fb_env, error)
    CALL fb_env_get(fb_env=fb_env, &
                    atomic_halos=atomic_halos, &
                    error=error)
    CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
                                 nhalos=nhalos, &
                                 halos=halos, &
                                 error=error)
    DO ispin = 1, nspin
       DO jhalo = 1, nhalos
          CALL fb_atomic_halo_get(atomic_halo=halos(jhalo), &
                                  owner_atom=jatom, &
                                  owner_id_in_halo=jatom_in_halo, &
                                  natoms=natoms_in_halo, &
                                  halo_atoms=halo_atoms, &
                                  errpr=error)
          ! construct atomic matrix for H for this_atom
          ALLOCATE(atomic_H_blk_row_start(natoms_in_halo + 1), &
                   atomic_H_blk_col_start(natoms_in_halo + 1), &
                   STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL fb_calc_atomic_matrix_size(mat_H(ispin)%matrix, &
                                          atomic_halo, &
                                          nrows, &
                                          ncols, &
                                          atomic_H_blk_row_start, &
                                          atomic_H_blk_col_start, &
                                          error)
          ALLOCATE(atomic_H(nrows,ncols), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL fb_build_atomic_matrix(mat_H(ispin)%matrix, &
                                      atomic_halo, &
                                      para_env, &
                                      atomic_H, &
                                      atomic_H_blk_row_start, &
                                      atomic_H_blk_col_start, &
                                      error)
          ! construct atomic matrix for S for this_atom
          ALLOCATE(atomic_S_blk_row_start(natoms_in_halo + 1), &
                   atomic_S_blk_col_start(natoms_in_halo + 1), &
                   STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL fb_calc_atomic_matrix_size(mat_S(ispin)%matrix, &
                                          atomic_halo, &
                                          nrows, &
                                          ncols, &
                                          atomic_S_blk_row_start, &
                                          atomic_S_blk_col_start, &
                                          error)
          ALLOCATE(atomic_S(nrows,ncols), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL fb_build_atomic_matrix(mat_S(ispin)%matrix, &
                                      atomic_halo, &
                                      para_env, &
                                      atomic_S, &
                                      atomic_S_blk_row_start, &
                                      atomic_S_blk_col_start, &
                                      error)

          ! diagonalise using LAPACK
          mat_dim = SIZE(atomic_H, 1)
          ALLOCATE(eigenvalues(mat_dim), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ALLOCATE(work(1), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ! get optimal work array size
          CALL dsygv(ITYPE = 1,           &
                     JOBZ  = 'V',         &
                     UPLO  = 'U',         &
                     N     = mat_dim,     &
                     A     = atomic_H,    &
                     LDA   = mat_dim,     &
                     B     = atomic_S,    &
                     LDB   = mat_dim,     &
                     W     = eigenvalues, &
                     WORK  = work,        &
                     LWORK = -1,          &
                     INFO  = info)
          work_array_size = work(1)
          ! now allocate work array
          DEALLOCATE(work, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ALLOCATE(work(work_array_size), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ! do calculation
          CALL dsygv(ITYPE = 1,               &
                     JOBZ  = 'V',             &
                     UPLO  = 'U',             &
                     N     = mat_dim,         &
                     A     = atomic_H,        &
                     LDA   = mat_dim,         &
                     B     = atomic_S,        &
                     LDB   = mat_dim,         &
                     W     = eigenvalues,     &
                     WORK  = work,            &
                     LWORK = work_array_size, &
                     INFO  = info)
          DEALLOCATE(work, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

          ! construct atomic filter matrix
          ! atomic_H now contains the eigenvectors
          
          ! first get the filter function.
          ALLOCATE(filter_function(mat_dim), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ! get system mu
          CALL get_mo_set(mo_set=mo_list(ispin), &
                          mu=fermi_level, &
                          maxocc=maxocc)
          CALL fb_filter_fermi_dirac_mu(filter_function, &
                                        eigenvalues, &
                                        filter_temperature, &
                                        fermi_level, &
                                        maxocc)
          ! atomic_H now has the eigenvectors, construct the
          ! version of it filtered through the filter function
          ALLOCATE(filtered_eigenvectors(mat_dim,mat_dim), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          DO jj = 1, mat_dim
             DO ii = 1, mat_dim
                filtered_eigenvectors(ii,jj) = &
                     filter_function(jj) * atomic_H(ii,jj)
             END DO ! ii
          END DO ! jj
          ! construct atomic filter matrix
          ALLOCATE(atomic_filter_matrix(mat_dim,mat_dim), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          ALLOCATE(temp_matrix(mat_dim,mat_dim), STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL dgemm(TRANSA = "N",                   &
                     TRANSB = "T",                   &
                     M      = mat_dim,               &
                     N      = mat_dim,               &
                     K      = mat_dim,               &
                     ALPHA  = 1.0_dp,                &
                     A      = filtered_eigenvectors, &
                     LDA    = mat_dim,               &
                     B      = atomic_H,              &
                     LDB    = mat_dim,               &
                     BETA   = 0.0_dp,                &
                     C      = atomic_filter_matrix,  &
                     LDC    = mat_dim)
          
          ! construct filter matrix
          
       


          ! now construct the local filter matrix blocks
          DO ihalo_atom = 1, natoms_in_halo
             iatom_global = halo_atoms(ihalo_atom)
             nrows_in_blk = atomic_H_blk_row_start(ihalo_atom + 1) - &
                            atomic_H_blk_row_start(ihalo_atom)
              DO jhalo_atom = 1, natoms_in_halo
                jatom_global = halo_atoms(jhalo_atom)
                atomic_kind => particle_set(jatom_global)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind, &
                                     kind_number=jkind)
                ncols_in_blk = nfns(jkind)
                ALLOCATE(filter_matrix_blk(nrows_in_blk,ncols_in_blk), STAT=stat)
                CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
                filter_matrix_blk = 0.0_dp
                DO itrial = 1, nfns(jkind)
                   CALL dgemv(                                                 &
                        TRANSA = "N",                                          &
                        TRANSB = "N",                                          &
                        M      = nrows_in_blk,                                 &
                        N      = mat_dim,                                      &
                        ALPHA  = 1.0_dp,                                       &
                        A      = atomic_filter_matrix(                         &
                                   atomic_H_blk_row_start(ihalo_atom) :        &
                                   atomic_H_blk_row_start(ihalo_atom + 1) - 1, &
                                   1 : mat_dim                                 &
                                 ),                                            &
                        LDA    = nrows_in_blk,                                 &
                        X      = atomic_S(                                     &
                                   1 : mat_dim,                                &
                                   atomic_S_blk_row_start(jhalo_atom) +        &
                                   fns(itrial,jkind) - 1                       &
                                 ),                                            &
                        INCX   = 1,                                            &
                        BETA   = 0.0_dp,                                       &
                        Y      = filter_matrix_blk(                            &
                                   1 : nrows_in_blk,                           &
                                   itrial                                      &
                                 ),                                            &
                        INCY   = 1)
                END DO ! itrial
                
                ! communicate the block
                
                ! construct communication tasks
                CALL fb_com_tasks_create(com_tasks_send, error)
                CALL fb_com_tasks_create(com_tasks_recv, error)
                stored_row = iatom_global
                stored_col = jatom_global
                ! number of send tasks will not exceed the total
                ! number of blocks generated on this process
                ntasks_send = 

                ! source is always this process
                src = para_env%mepos
                ! destination is the proc this block belongs to
                CALL cp_dbcsr_get_stored_coordinates(matrix=filter_matrix(ispin)%matrix, &
                                                     row=stored_row, col=stored_col, &
                                                     transpose=.FALSE., &
                                                     processor=dest)
                


                ! put the filter matrix block into the filter matrix
                CALL cp_dbcsr_put_block(matrix=filter_matrix(ispin)%matrix, &
                                        row=iatom_global, &
                                        col=jatom_global, &
                                        block=filter_matrix_blk)
                

             END DO ! jhalo_atoms
          END DO ! ihalo_atom

       END DO  ! ihalo
    END DO ! ispin

  END SUBROUTINE fb_env_build_filter_matrix


  SUBROUTINE fb_filter_fermi_dirac_mu(f, eigenvals, T, mu, maxocc)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in) :: T
    REAL(KIND=dp), INTENT(in) :: mu
    REAL(KIND=dp), INTENT(in) :: maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_filter_fermi_dirac_mu', &
                                   routineP = moduleN//':'//routineN
    
    LOGICAL :: failure
    INTEGER :: ne
    REAL(KIND=dp) :: kTS

    CALL Fermi(f, ne, kTS, eigenvals, mu, T, maxocc)    
  END SUBROUTINE fb_filter_fermi_dirac_mu


  SUBROUTINE fb_filter_fermi_dirac_ne(f, eigenvals, T, ne, maxocc)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in) :: T
    REAL(KIND=dp), INTENT(in) :: ne
    REAL(KIND=dp), INTENT(in) :: maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_filter_fermi_dirac_ne', &
                                   routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: mu, kTS
    ! mu is the calculated fermi level
    ! kTS is the calculated entropic contribution to the energy i.e. -TS
    ! kTS= kT*[f ln f + (1-f) ln (1-f)]
    LOGICAL       :: failure

    failure = .FALSE.    
    CALL FermiFixed(f, mu, kTS, eigenval, ne, fermi_T, maxocc)
  END SUBROUTINE fb_filter_fermi_dirac_ne


END MODULE qs_fb_env_methods
