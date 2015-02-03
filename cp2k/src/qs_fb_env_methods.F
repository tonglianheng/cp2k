MODULE qs_fb_env_methods

  USE kinds, ONLY: dp, &
                   default_string_length
  USE qs_fb_env_types, ONLY: fb_env_obj, &
                             fb_env_get, &
                             fb_env_set, &
                             fb_env_has_data
  USE qs_environment_types,    ONLY: qs_environment_type, &
                                     get_qs_env
!BEG:DEBUG:LT:2015/02/02
!  USE qs_environment_methods, ONLY: qs_create_task_list
!END:DEBUG:LT:2015/02/02
  USE qs_fb_atomic_halo_types, ONLY: fb_atomic_halo_obj, &
                                     fb_atomic_halo_nullify, &
                                     fb_atomic_halo_create, &
                                     fb_atomic_halo_get, &
                                     fb_atomic_halo_set, &
                                     fb_atomic_halo_sort, &
                                     fb_atomic_halo_nelectrons_estimate_Z, &
                                     fb_atomic_halo_list_obj, &
                                     fb_atomic_halo_list_init, &
                                     fb_atomic_halo_list_create, &
                                     fb_atomic_halo_list_get, &
                                     fb_atomic_halo_list_set, &
                                     fb_atomic_halo_list_has_data, &
                                     fb_atomic_halo_list_nullify
  USE particle_types,        ONLY: particle_type
  USE atomic_kind_types,     ONLY: atomic_kind_type, &
                                   get_atomic_kind
  USE qs_neighbor_lists,      ONLY: pair_radius_setup
  USE distribution_1d_types, ONLY: distribution_1d_type
  USE cp_dbcsr_interface, ONLY: cp_dbcsr_type, &
                                cp_dbcsr_p_type, &
                                cp_dbcsr_allocate_matrix_set, &
                                cp_dbcsr_init, &
                                cp_dbcsr_create, &
                                cp_dbcsr_distribution, &
                                dbcsr_type_no_symmetry, &
                                cp_dbcsr_row_block_sizes, &
                                cp_dbcsr_col_block_sizes, &
                                cp_dbcsr_nfullrows_total, &
                                cp_dbcsr_nfullcols_total, &
                                cp_dbcsr_finalize, &
                                cp_dbcsr_multiply, &
                                cp_dbcsr_release, &
                                cp_dbcsr_set, &
                                cp_dbcsr_reserve_blocks, &
                                cp_dbcsr_iterator, &
                                cp_dbcsr_iterator_start, &
                                cp_dbcsr_iterator_stop, &
                                cp_dbcsr_iterator_next_block, &
                                cp_dbcsr_iterator_blocks_left, &
                                cp_dbcsr_get_matrix_type, &
                                array_i1d_obj, &
!BEG:DEBUG:LT:2015/01/20
                                cp_dbcsr_print, &
                                cp_dbcsr_desymmetrize, &
!END:DEBUG:LT:2015/01/20
                                array_nullify
  USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                 cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
  USE input_constants, ONLY: use_aux_fit_basis_set,&
                             use_orb_basis_set, &
                             cholesky_inverse
  USE orbital_pointers, ONLY: ncoset, &
                              nco
  USE qs_mo_types, ONLY: mo_set_p_type, &
                         get_mo_set, &
                         set_mo_set, &
                         init_mo_set, &
                         set_mo_occupation, &
                         allocate_mo_set, &
                         deallocate_mo_set
  USE qs_mo_methods, ONLY: calculate_density_matrix
  USE string_utilities, ONLY: compress,&
                              uppercase
  USE qs_diis, ONLY: qs_diis_b_step
  USE qs_scf_types, ONLY: qs_scf_env_type
  USE scf_control_types, ONLY: scf_control_type
  USE input_section_types, ONLY: section_vals_type, &
                                 section_vals_val_get, &
                                 section_vals_get_subs_vals
  USE cp_para_types, ONLY: cp_para_env_type
  USE cp_blacs_env, ONLY: cp_blacs_env_type
  USE cp_fm_types, ONLY: cp_fm_type, &
                         cp_fm_create, &
                         cp_fm_set_all, &
                         cp_fm_release, &
                         cp_fm_to_fm, &
!BEG:DEBUG:LT:2015/02/02
                         cp_fm_write
!END:DEBUG:LT:2015/02/02
  USE qs_fb_trial_fns_types, ONLY: fb_trial_fns_obj, &
                                   fb_trial_fns_nullify, &
                                   fb_trial_fns_get, &
                                   fb_trial_fns_set, &
                                   fb_trial_fns_create
  USE qs_fb_filter_matrix_methods, ONLY: fb_fltrmat_build
  USE qs_matrix_pools, ONLY: qs_matrix_pools_type, &
                             mpools_create, &
                             mpools_rebuild_fm_pools, &
                             mpools_release
  USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                          cp_fm_struct_type, &
                          cp_fm_struct_release
  USE cp_fm_cholesky, ONLY: cp_fm_cholesky_decompose
  USE cp_fm_basic_linalg, ONLY: cp_fm_triangular_invert
  USE qs_scf_methods, ONLY: eigensolver
  USE basis_set_types, ONLY: gto_basis_set_type, &
                             gto_basis_set_p_type, &
                             get_gto_basis_set
  USE qs_integral_utils, ONLY: basis_set_list_setup
  USE qs_kind_types, ONLY: qs_kind_type, &
                           get_qs_kind
  USE message_passing, ONLY: mp_max
  USE cp_fm_basic_linalg, ONLY: cp_fm_gemm
  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type


#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_methods'

  PUBLIC :: fb_env_do_diag, &
            fb_env_read_input, &
            fb_env_build_rcut_auto, &
            fb_env_build_atomic_halos, &
            fb_env_build_trial_fns_auto


CONTAINS


  SUBROUTINE fb_env_do_diag(fb_env, qs_env, matrix_ks, matrix_s, scf_section, &
                            diis_step, error)
    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_s
    TYPE(section_vals_type), POINTER :: scf_section
    LOGICAL, INTENT(INOUT) :: diis_step
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_env_do_diag', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure, use_jacobi
    INTEGER :: ispin, nspin, stat, &
               original_nfullrowsORcols_total, &
               filtered_nfullrowsORcols_total, &
               nelectron, nfullrows_total, nfullcols_total, &
               homo_filtered, lfomo_filtered, &
               nmo_filtered, max_nelectron
    REAL(kind=dp) :: diis_error, eps_diis
    REAL(kind=dp) :: fermi_level, maxocc, flexible_electron_count, &
                     filter_temp, mu_filtered, KTS_filtered
    REAL(KIND=dp), DIMENSION(:), POINTER :: eigenvalues_filtered, occ_filtered, &
                                            eigenvalues, occ
    TYPE(cp_fm_type), POINTER :: mo_coeff_filtered, mo_coeff, &
                                 fm_matrix_filtered_s, fm_matrix_filtered_ks, &
                                 fm_matrix_work, fm_matrix_filter
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos, mos_filtered
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_filter
    TYPE(cp_dbcsr_type) :: matrix_tmp, matrix_filtered_ks, &
                           matrix_filtered_s
!BEG:DEBUG:LT:2015/01/20
    TYPE(cp_dbcsr_p_type) :: matrix_ks_desymm, matrix_s_desymm
!END:DEBUG:LT:2015/01/20
    TYPE(cp_dbcsr_type), POINTER :: matrix_filtered_p
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_trial_fns_obj) :: trial_fns
    CHARACTER(len=2) :: spin_string
    CHARACTER(len=default_string_length) :: name
    TYPE(array_i1d_obj) :: filtered_rowORcol_block_sizes, &
                           original_rowORcol_block_sizes
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(qs_matrix_pools_type), POINTER :: my_mpools
    TYPE(cp_fm_struct_type), POINTER  :: fm_struct
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb

!BEG:DEBUG:LT:2015/01/29
    INTEGER :: LT_UNIT
!END:DEBUG:LT:2015/01/29

    failure = .FALSE.

    NULLIFY(scf_env, scf_control, para_env, blacs_env, mos, sab_orb, &
            mos_filtered, matrix_filter, particle_set, &
            matrix_filtered_p, eigenvalues_filtered, occ_filtered, &
            mo_coeff_filtered, eigenvalues, occ, mo_coeff)
    CALL fb_atomic_halo_list_nullify(atomic_halos)
    CALL fb_trial_fns_nullify(trial_fns)
    CALL array_nullify(original_rowORcol_block_sizes)
    CALL array_nullify(filtered_rowORcol_block_sizes)

    nspin = SIZE(matrix_ks)

    ! get qs_env information
    CALL get_qs_env(qs_env=qs_env, &
                    scf_env=scf_env, &
                    scf_control=scf_control, &
                    para_env=para_env, &
                    blacs_env=blacs_env, &
                    particle_set=particle_set, &
                    mos=mos, &
                    error=error)

!BEG:DEBUG:LT:2015/01/30
    LT_UNIT = 100 + para_env%mepos
!END:DEBUG:LT:2015/01/30

    DO ispin = 1, nspin
       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,&
                             scf_env%scf_work1(ispin)%matrix,&
                             error=error)
    END DO
    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
       CALL qs_diis_b_step(scf_env%scf_diis_buffer, mos, scf_env%scf_work1, &
                           scf_env%scf_work2, scf_env%iter_delta, &
                           diis_error, diis_step, eps_diis, scf_control%nmixing, &
                           s_matrix=matrix_s, scf_section=scf_section, error=error)
    ELSE
       diis_step = .FALSE.
    END IF

    IF ((scf_env%iter_count > 1).AND.&
        (scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    IF (diis_step) THEN
       scf_env%iter_param = diis_error
       IF (use_jacobi) THEN
          scf_env%iter_method = "DIIS/Jacobi"
       ELSE
          scf_env%iter_method = "DIIS/Diag."
       END IF
    ELSE
       IF(scf_env%mixing_method ==0) THEN
          scf_env%iter_method = "NoMix/Diag."
       ELSE IF(scf_env%mixing_method ==1) THEN
          scf_env%iter_param = scf_env%p_mix_alpha
          IF (use_jacobi) THEN
             scf_env%iter_method = "P_Mix/Jacobi"
          ELSE
             scf_env%iter_method = "P_Mix/Diag."
          END IF
       ELSEIF(scf_env%mixing_method>1) THEN
          scf_env%iter_param = scf_env%mixing_store%alpha
          IF (use_jacobi) THEN
             scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Jacobi"
          ELSE
             scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
          END IF
       END IF
    END IF


    ! get filter temperature
    CALL fb_env_get(fb_env=fb_env, &
                    filter_temperature=filter_temp, &
                    error=error)

    ! construct atomic halos: atomic halos should be constructed
    ! here because atoms would have moved at every MD step.
    CALL fb_env_build_atomic_halos(fb_env, qs_env, error)
    CALL fb_env_get(fb_env=fb_env, &
                    atomic_halos=atomic_halos, &
                    error=error)
!TODO:LT:2015/01/19 can move the above into the outer MD loop sections?

    ! construct trial functions
    max_nelectron = 0
    DO ispin = 1, nspin
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       nelectron=nelectron, &
                       maxocc=maxocc)
       max_nelectron = MAX(max_nelectron, nelectron)
    END DO
    CALL fb_env_build_trial_fns_auto(fb_env, qs_env, REAL(max_nelectron,dp), &
                                     maxocc, error)
    CALL fb_env_get(fb_env=fb_env, &
                    trial_fns=trial_fns, &
                    error=error)

    ! allocate filter matrix (matrix_filter(ispin)%matrix are
    ! nullified by cp_dbcsr_allocate_matrix_set)
    CALL cp_dbcsr_allocate_matrix_set(matrix_filter, nspin, error)
    DO ispin = 1, nspin
       ! get system-wide fermi energy and occupancy, we use this to
       ! define the filter function used for the filter matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       mu=fermi_level, &
                       maxocc=maxocc)

!BEG:DEBUG:LT:2015/01/2

       ! desymmetrise matrix first to avoid the complex and hard to
       ! understand DBCSR symmetric layout first.
       NULLIFY(matrix_ks_desymm%matrix)
       ALLOCATE(matrix_ks_desymm%matrix)
       CALL cp_dbcsr_init(matrix_ks_desymm%matrix, error)
       CALL cp_dbcsr_desymmetrize(matrix_ks(ispin)%matrix, matrix_ks_desymm%matrix, error)
       NULLIFY(matrix_s_desymm%matrix)
       ALLOCATE(matrix_s_desymm%matrix)
       CALL cp_dbcsr_init(matrix_s_desymm%matrix, error)
       CALL cp_dbcsr_desymmetrize(matrix_s(ispin)%matrix, matrix_s_desymm%matrix, error)

!END:DEBUG:LT:2015/01/2

       ! get filter matrix name
       WRITE (spin_string, FMT="(I1)"), ispin
       name = TRIM("FILTER MATRIX SPIN "//spin_string)
       CALL compress(name)
       CALL uppercase(name)
       ! calculate filter matrix (matrix_s(1) is the overlap, the rest
       ! in the array are its derivatives)
!BEG:ORIG:LT:2015/01/2
       ! CALL fb_fltrmat_build(H_mat=matrix_ks(ispin)%matrix, &
       !                       S_mat=matrix_s(1)%matrix, &
       !                       atomic_halos=atomic_halos, &
       !                       trial_fns=trial_fns, &
       !                       para_env=para_env, &
       !                       particle_set=particle_set, &
       !                       fermi_level=fermi_level, &
       !                       max_occ=maxocc, &
       !                       filter_temp=filter_temp, &
       !                       name=name, &
       !                       filter_mat=matrix_filter(ispin)%matrix, &
       !                       error=error)
!END:ORIG:LT:2015/01/2

!BEG:DEBUG:LT:2015/01/30
       write (LT_UNIT,*) "|LT| fb_env_do_diag: atomic_halos associated = ", &
             fb_atomic_halo_list_has_data(atomic_halos)
!END:DEBUG:LT:2015/01/30


       CALL fb_fltrmat_build(H_mat=matrix_ks_desymm%matrix, &
                             S_mat=matrix_s_desymm%matrix, &
                             atomic_halos=atomic_halos, &
                             trial_fns=trial_fns, &
                             para_env=para_env, &
                             particle_set=particle_set, &
                             fermi_level=fermi_level, &
                             max_occ=maxocc, &
                             filter_temp=filter_temp, &
                             name=name, &
                             filter_mat=matrix_filter(ispin)%matrix, &
                             error=error)
!BEG:DEBUG:LT:2015/01/2
       CALL cp_dbcsr_release(matrix_ks_desymm%matrix, error)
       CALL cp_dbcsr_release(matrix_s_desymm%matrix, error)
       DEALLOCATE(matrix_ks_desymm%matrix)
       DEALLOCATE(matrix_s_desymm%matrix)
!END:DEBUG:LT:2015/01/2

    END DO

    ! Obtain matrix dimensions. KS and S matrices are symmetric, so
    ! row_block_sizes and col_block_sizes should be identical. The
    ! same applies to the filtered block sizes. Note that filter
    ! matrix will have row_block_sizes equal to that of the original,
    ! and col_block_sizes equal to that of the filtered.  We assume
    ! also that the matrix dimensions are identical for both spin
    ! channels.
    original_rowORcol_block_sizes = cp_dbcsr_row_block_sizes(matrix_ks(1)%matrix)
    filtered_rowORcol_block_sizes = cp_dbcsr_col_block_sizes(matrix_filter(1)%matrix)
    original_nfullrowsORcols_total = cp_dbcsr_nfullrows_total(matrix_ks(1)%matrix)
    filtered_nfullrowsORcols_total = cp_dbcsr_nfullcols_total(matrix_filter(1)%matrix)

    ! filter diagonalisation works on a smaller basis set, and thus
    ! requires a new mo_set (molecular orbitals | eigenvectors) and
    ! the corresponding matrix pools for the eigenvector coefficients
    ALLOCATE(mos_filtered(nspin), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ispin = 1, nspin
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       maxocc=maxocc, &
                       nelectron=nelectron, &
                       flexible_electron_count=flexible_electron_count)
       NULLIFY(mos_filtered(ispin)%mo_set)
       CALL allocate_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                            nao=filtered_nfullrowsORcols_total, &
                            nmo=filtered_nfullrowsORcols_total, &
                            nelectron=nelectron, &
                            n_el_f=REAL(nelectron, dp), &
                            maxocc=maxocc, &
                            flexible_electron_count=flexible_electron_count, &
                            error=error)
    END DO
    NULLIFY(my_mpools)
    CALL mpools_create(mpools=my_mpools, error=error)
    CALL mpools_rebuild_fm_pools(mpools=my_mpools,&
                                 mos=mos_filtered,&
                                 blacs_env=blacs_env,&
                                 para_env=para_env,&
                                 error=error)

    ! create DBCSR filtered KS matrix, this is reused for each spin
    ! channel
    CALL cp_dbcsr_init(matrix_filtered_ks, error=error)
    ! both row_blk_size and col_blk_size should be that of
    ! col_blk_size of the filter matrix
    CALL cp_dbcsr_create(matrix=matrix_filtered_ks, &
                         name=TRIM("FILTERED_KS_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_ks(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_ks, error=error)

    ! create DBCSR filtered S (overlap) matrix. Note that
    ! matrix_s(1)%matrix is the orginal overlap matrix---the rest in
    ! the array are derivatives, and it should not depend on
    ! spin. HOWEVER, since the filter matrix is constructed from KS
    ! matrix, and does depend on spin, the filtered S also becomes
    ! spin dependent. Nevertheless this matrix is reused for each spin
    ! channel
    CALL cp_dbcsr_init(matrix_filtered_s, error=error)
    ! both row_blk_size and col_blk_size should be that of
    ! col_blk_size of the filter matrix
    CALL cp_dbcsr_create(matrix=matrix_filtered_s, &
                         name=TRIM("FILTERED_S_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_s(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_s, error=error)

    ! create temporary matrix for constructing filtered KS and S
    CALL cp_dbcsr_init(matrix_tmp, error=error)
    ! the temporary matrix won't be square
    CALL cp_dbcsr_create(matrix=matrix_tmp, &
                         name=TRIM("TEMPORARY_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_s(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=original_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_tmp, error=error)

    ! create fm format matrices used for diagonalisation
!BEG:DEBUG:LT:2015/01/20
    NULLIFY(fm_matrix_filter)
!END:DEBUG:LT:2015/01/20
    NULLIFY(fm_matrix_filtered_s)
    NULLIFY(fm_matrix_filtered_ks)
    NULLIFY(fm_matrix_work)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fmstruct=fm_struct, &
                             para_env=para_env, &
                             context=blacs_env, &
                             nrow_global=filtered_nfullrowsORcols_total, &
                             ncol_global=filtered_nfullrowsORcols_total, &
                             error=error)
    ! both fm_matrix_filtered_s and fm_matrix_filtered_ks are reused
    ! for each spin channel
!BEG:DEBUG:LT:2015/01/20
    CALL cp_fm_create(fm_matrix_filter, &
                      fm_struct, &
                      name="FM_MATRIX_FILTER", &
                      error=error)
!END:DEBUG:LT:2015/01/20
    CALL cp_fm_create(fm_matrix_filtered_s, &
                      fm_struct, &
                      name="FM_MATRIX_FILTERED_S", &
                      error=error)
    CALL cp_fm_create(fm_matrix_filtered_ks, &
                      fm_struct, &
                      name="FM_MATRIX_FILTERED_KS", &
                      error=error)
    ! creaate work matrix
    CALL cp_fm_create(fm_matrix_work, fm_struct, name="FM_MATRIX_WORK", error=error)
    ! all fm matrices are created, so can release fm_struct
    CALL cp_fm_struct_release(fm_struct, error=error)

    ! construct filtered KS, S matrix and diagonalise
    DO ispin = 1, nspin

       ! construct filtered KS matrix
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_ks(ispin)%matrix, matrix_filter(ispin)%matrix, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_tmp, &
                              0.0_dp, matrix_filtered_ks, error=error)
       ! construct filtered S_matrix
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_s(1)%matrix, matrix_filter(ispin)%matrix, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_tmp, &
                              0.0_dp, matrix_filtered_s, error=error)

       ! now that we have the filtered KS and S matrices for this spin
       ! channel, perform ordinary diagonalisation

       ! convert DBCSR matrices to fm format
       CALL copy_dbcsr_to_fm(matrix_filtered_s, fm_matrix_filtered_s, error=error)
       CALL copy_dbcsr_to_fm(matrix_filtered_ks, fm_matrix_filtered_ks, error=error)

!BEG:DEBUG:LT:2015/01/20
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_filter(1)%matrix: (DBCSR)"
       CALL cp_dbcsr_print(matrix_filter(1)%matrix, unit_nr=LT_UNIT, error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_s(1)%matrix: (DBCSR):"
       CALL cp_dbcsr_print(matrix_s(1)%matrix, unit_nr=LT_UNIT, error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_filtered_s(1)%matrix: (DBCSR)"
       CALL cp_dbcsr_print(matrix_filtered_s, unit_nr=LT_UNIT, error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: fm_matrix_filtered_s: (FM)"
       CALL cp_fm_write(matrix=fm_matrix_filtered_s, &
                        long_description=.TRUE., &
                        local=.TRUE., &
                        unit_nr=LT_UNIT, &
                        error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_ks(1)%matrix: (DBCSR):"
       CALL cp_dbcsr_print(matrix_ks(1)%matrix, unit_nr=LT_UNIT, error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_filtered_ks(1)%matrix: (DBCSR)"
       CALL cp_dbcsr_print(matrix_filtered_ks, unit_nr=LT_UNIT, error=error)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: fm_matrix_filtered_ks: (FM)"
       CALL cp_fm_write(matrix=fm_matrix_filtered_ks, &
                        long_description=.TRUE., &
                        local=.TRUE., &
                        unit_nr=LT_UNIT, &
                        error=error)
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/01/20

       ! now diagonalise

       ! calculate S^(-1/2) (cholescky decomposition)
       CALL cp_fm_cholesky_decompose(fm_matrix_filtered_s, error=error)
       CALL cp_fm_triangular_invert(fm_matrix_filtered_s, error=error)
       ! setup matrix pools for the molecular orbitals
       CALL init_mo_set(mos_filtered(ispin)%mo_set,&
                        my_mpools%ao_mo_fm_pools(ispin)%pool,&
                        name="FILTERED_MOS",&
                        error=error)
       ! initialise work matrix
       CALL cp_fm_set_all(matrix=fm_matrix_work, alpha=0.0_dp, error=error)
       ! diagonalise

!BEG:DEBUG:LT:2015/01/10
       scf_env%cholesky_method = cholesky_inverse
       use_jacobi = .FALSE.
!END:DEBUG:LT:2015/01/10

       CALL eigensolver(matrix_ks_fm=fm_matrix_filtered_ks, &
                        mo_set=mos_filtered(ispin)%mo_set, &
                        ortho=fm_matrix_filtered_s, &
                        work=fm_matrix_work,&
                        cholesky_method=scf_env%cholesky_method, &
                        use_jacobi=use_jacobi, &
                        error=error)

!BEG:DEBUG:LT:2015/02/02
       CALL get_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                       eigenvalues=eigenvalues_filtered, &
                       mo_coeff=mo_coeff_filtered)
       write (LT_UNIT,*) "|LT| fb_env_do_diag: eigenvalues_filtered = "
       write (LT_UNIT,*) eigenvalues_filtered
       write (LT_UNIT,*) "|LT| fb_env_do_diag: mo_coeff_filtered: (FM)"
       CALL cp_fm_write(matrix=mo_coeff_filtered, &
                        long_description=.TRUE., &
                        local=.TRUE., &
                        unit_nr=LT_UNIT, &
                        error=error)
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02

    END DO ! spin

    ! release temporary matrices
    CALL cp_dbcsr_release(matrix_filtered_s, error)
    CALL cp_dbcsr_release(matrix_filtered_ks, error)
    CALL cp_fm_release(fm_matrix_filtered_s, error=error)
    CALL cp_fm_release(fm_matrix_filtered_ks, error=error)
    CALL cp_fm_release(fm_matrix_work, error=error)

    ! calculate filtered molecular orbital occupation numbers and fermi
    ! level etc
    CALL set_mo_occupation(mo_array=mos_filtered, &
                           smear=scf_control%smear, &
                           error=error)

!BEG:DEBUG:LT:2015/02/02
    CALL get_mo_set(mo_set=mos_filtered(1)%mo_set, &
                    occupation_numbers=occ_filtered, &
                    mu=mu_filtered)
    write (LT_UNIT,*) "|LT| fb_env_do_diag: associated(occ_filtered) = ", associated(occ_filtered)
    CALL FLUSH(LT_UNIT)
    write (LT_UNIT,*) "|LT| fb_env_do_diag: SIZE(occ_filtered) = ", SIZE(occ_filtered)
    CALL FLUSH(LT_UNIT)
    write (LT_UNIT,*) "|LT| fb_env_do_diag: occ_filtered = ", occ_filtered
    write (LT_UNIT,*) "|LT| fb_env_do_diag: mu_filtered = ", mu_filtered
    CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02


    ! get the filtered density matrix and then convert back to the
    ! full basis version in scf_env ready to be used outside this
    ! subroutine
    ALLOCATE(matrix_filtered_p, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL cp_dbcsr_init(matrix_filtered_p, error=error)
    ! the filtered density matrix should have the same sparse
    ! structure as the original density matrix, we must copy the
    ! sparse structure here, since construction of the density matrix
    ! preserves its sparse form, and therefore matrix_filtered_p must
    ! have its blocks allocated here now. We assume the original
    ! density matrix scf_env%p_mix_new has the same sparse structure
    ! in both spin channels.
    CALL cp_dbcsr_create(matrix=matrix_filtered_p, &
                         name=TRIM("FILTERED_MATRIX_P"), &
                         dist=cp_dbcsr_distribution(scf_env%p_mix_new(1)%matrix), &
                         matrix_type=cp_dbcsr_get_matrix_type(scf_env%p_mix_new(1)%matrix), &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_p, error=error)
    CALL fb_dbcsr_copy_sparse_struct(matrix_filtered_p, &
                                     scf_env%p_mix_new(1)%matrix, &
                                     error)
    ! old implementation, using sab_orb to allocate the blocks in matrix_filtered_p
    ! CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, error=error)
    ! CALL cp_dbcsr_alloc_block_from_nbl(matrix_filtered_p, sab_orb, error=error)
    CALL cp_dbcsr_set(matrix_filtered_p, 0.0_dp, error=error)

    DO ispin = 1, nspin

!BEG:DEBUG:LT:2015/02/02
       write (LT_UNIT,*) "|LT| fb_env_do_diag: before update density: scf_env%p_mix_new: (DBCSR)"
       CALL cp_dbcsr_print(scf_env%p_mix_new(ispin)%matrix, unit_nr=LT_UNIT, error=error)
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02

      ! calculate matrix_filtered_p
       CALL calculate_density_matrix(mos_filtered(ispin)%mo_set, &
                                     matrix_filtered_p, &
                                     error=error)


!BEG:DEBUG:LT:2015/02/02
       write (LT_UNIT,*) "|LT| fb_env_do_diag: matrix_filtered_p: (DBCSR) "
       CALL cp_dbcsr_print(matrix_filtered_p, unit_nr=LT_UNIT, error=error)
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02


       ! convert back to full basis p
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_filtered_p, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                              matrix_tmp, matrix_filter(ispin)%matrix, &
                              0.0_dp, scf_env%p_mix_new(ispin)%matrix, error=error)

!BEG:DEBUG:LT:2015/02/02
       write (LT_UNIT,*) "|LT| fb_env_do_diag: after update density: scf_env%p_mix_new: (DBCSR)"
       CALL cp_dbcsr_print(scf_env%p_mix_new(ispin)%matrix, unit_nr=LT_UNIT, error=error)
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02


    END DO

    ! release temporary matrices
    CALL cp_dbcsr_release(matrix_tmp, error)
    CALL cp_dbcsr_release(matrix_filtered_p, error)
    DEALLOCATE(matrix_filtered_p, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! we still need to convert mos_filtered back to the full basis
    ! version (mos) for this, we need to update mo_coeff (and/or
    ! mo_coeff_b --- the DBCSR version, if used) of mos

    ! note also that mo_eigenvalues cannot be fully updated, given
    ! that the eigenvalues are computed in a smaller basis, and thus
    ! do not give the full spectron. Printing of molecular states
    ! (molecular DOS) at each SCF step is therefore not recommended
    ! when using this method. The idea is that if one wants a full
    ! molecular DOS, then one should perform a full diagonalisation
    ! without the filters once the SCF has been achieved.

    ! NOTE: from reading the source code, it appears that mo_coeff_b
    ! is actually never used by default (DOUBLE CHECK?!). Even
    ! subroutine eigensolver_dbcsr updates mo_coeff, and not
    ! mo_coeff_b.
    DO ispin = 1, nspin
       ! now the full basis mo_set should only contain the reduced
       ! number of eigenvectors and eigenvalues
       CALL get_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                       homo=homo_filtered, &
                       lfomo=lfomo_filtered, &
                       nmo=nmo_filtered, &
                       eigenvalues=eigenvalues_filtered, &
                       occupation_numbers=occ_filtered, &
                       mo_coeff=mo_coeff_filtered, &
                       kTS=kTS_filtered, &
                       mu=mu_filtered)
       CALL set_mo_set(mo_set=mos(ispin)%mo_set, &
                       homo=homo_filtered, &
                       lfomo=lfomo_filtered, &
                       nmo=nmo_filtered, &
                       kTS=kTS_filtered, &
                       mu=mu_filtered, &
                       error=error)
       ! CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
       !                 occupation_numbers=occ, &
       !                 eigenvalues=eigenvalues, &
       !                 mo_coeff=mo_coeff)
       ! eigenvalues(:) = 0.0_dp
       ! eigenvalues(1:nmo_filtered) = eigenvalues_filtered(1:nmo_filtered)
       ! occ(:) = 0.0_dp
       ! occ(1:nmo_filtered) = occ_filtered(1:nmo_filtered)

!BEG:ORIG:LT:2015/01/20
       ! CALL cp_dbcsr_sm_fm_multiply(matrix_filter(ispin)%matrix, &
       !                              mo_coeff_filtered, &
       !                              mo_coeff, &
       !                              nmo_filtered, &
       !                              error=error)
!END:ORIG:LT:2015/01/20

!BEG:DEBUG:LT:2015/01/20
       ! NULLIFY(fm_matrix_work)
       ! CALL cp_fm_create(fm_matrix_work, &
       !                   mos_filtered(ispin)%mo_set%mo_coeff%matrix_struct, error=error)
       ! CALL cp_fm_to_fm(mos_filtered(ispin)%mo_set%mo_coeff,fm_matrix_work,error=error)
       ! CALL copy_dbcsr_to_fm(matrix_filter(ispin)%matrix, fm_matrix_filter, error=error)
       ! CALL cp_fm_gemm("N","N", &
       !                 original_nfullrowsORcols_total, &
       !                 filtered_nfullrowsORcols_total, &
       !                 filtered_nfullrowsORcols_total, &
       !                 1.0_dp, fm_matrix_filter, fm_matrix_work, &
       !                 0.0_dp, mo_coeff, error=error)
       ! CALL cp_fm_release(fm_matrix_work, error=error)
       ! CALL cp_fm_release(fm_matrix_filter, error=error)
!END:DEBUG:LT:2015/01/20
    END DO

    ! release temporary matrices
    CALL mpools_release(mpools=my_mpools, error=error)
    DO ispin = 1, nspin
       CALL deallocate_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                              error=error)
       CALL cp_dbcsr_release(matrix_filter(ispin)%matrix, error)
    END DO
    DEALLOCATE(mos_filtered, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(matrix_filter, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)


!BEG:DEBUG:LT:2015/02/02
!    CALL qs_create_task_list(qs_env, error)
!END:DEBUG:LT:2015/02/02

  END SUBROUTINE fb_env_do_diag


  SUBROUTINE fb_env_read_input(fb_env, scf_section, error)

    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(section_vals_type), POINTER :: scf_section
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_read_input', &
                                   routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER :: fb_section
    REAL(KIND=dp) :: filter_temperature

    NULLIFY(fb_section)

    fb_section => section_vals_get_subs_vals(scf_section, &
                                             "DIAGONALIZATION%FILTER_MATRIX", &
                                             error=error)
    CALL section_vals_val_get(fb_section, "FILTER_TEMPERATURE", &
                              r_val=filter_temperature, &
                              error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    filter_temperature=filter_temperature, &
                    error=error)

  END SUBROUTINE fb_env_read_input



  SUBROUTINE fb_env_build_rcut_auto(fb_env, qs_env, error)
    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_rcut_auto', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(qs_kind_type), DIMENSION(:), POINTER :: qs_kind_set
    TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: basis_set_list
    TYPE(gto_basis_set_type), POINTER :: basis_set

    REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
    REAL(KIND=dp) :: kind_radius
    INTEGER :: ikind, basis_set_id, nkinds, stat

    NULLIFY(rcut, qs_kind_set)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set, &
                    error=error)

    nkinds = SIZE(qs_kind_set)
    ALLOCATE(rcut(nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! reading from the other parts of the code, it seemed that
    ! aux_fit_basis_set is only used when do_admm is TRUE. This can be
    ! seen from the calls to generate_qs_task_list subroutine in
    ! qs_create_task_list, found in qs_environment_methods.F:
    ! basis_set_id is only set as input parameter for do_admm
    ! calculations, and if not set, the task list is generated using
    ! the default basis_set_id=use_orb_basis_set.
    IF (qs_env%dft_control%do_admm) THEN
       basis_set_id = use_aux_fit_basis_set
    ELSE
       basis_set_id = use_orb_basis_set
    END IF

    ALLOCATE(basis_set_list(nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL basis_set_list_setup(basis_set_list, basis_set_id, qs_kind_set)

    DO ikind = 1, nkinds
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              kind_radius=kind_radius)
       rcut(ikind) = 0.2_dp * kind_radius
    END DO

    CALL fb_env_set(fb_env=fb_env, &
                    rcut=rcut, &
                    error=error)

    ! cleanup
    DEALLOCATE(basis_set_list, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE fb_env_build_rcut_auto


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
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
    INTEGER :: stat, natoms_local, natoms_global, nkinds, ikind, &
               iatom, nkinds_global, jkind, jatom, iatom_local, &
               n_local_atoms_in_kind, n_global_atoms_in_kind, &
               jatom_in_kind, ii, ihalo, ihalo_atom, nhalo_atoms, &
               nelectrons_in_halo, owner_id_in_halo, max_natoms_local
    INTEGER, DIMENSION(:), POINTER :: halo_atoms
    INTEGER, DIMENSION(:), POINTER :: atom_list
    REAL(KIND=dp) :: rij
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: pair_radius
    LOGICAL, DIMENSION(:), ALLOCATABLE :: present_ab

    failure = .FALSE.

    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)

    NULLIFY(rcut, particle_set, atomic_kind_set, &
            local_particles, halos, halo_atoms, para_env)
    CALL fb_atomic_halo_list_nullify(atomic_halos)

    ! get relevant data from fb_env
    CALL fb_env_get(fb_env=fb_env, &
                    rcut=rcut, &
                    error=error)

    ! create atomic_halos
    CALL fb_atomic_halo_list_create(atomic_halos, error)

    ! get the number of atoms and kinds:
    CALL get_qs_env(qs_env=qs_env, &
                    local_particles=local_particles, &
                    natom=natoms_global, &
                    particle_set=particle_set, &
                    atomic_kind_set=atomic_kind_set, &
                    nkind=nkinds_global, &
                    para_env=para_env, &
                    error=error)

    nkinds = SIZE(local_particles%list)
    natoms_local = SUM(local_particles%n_el)
    ! get the maximum number of local atoms across the procs.
    max_natoms_local = natoms_local
    CALL mp_max(max_natoms_local, para_env%group)

!BEG:DEBUG:LT:2014/08/22
    CALL cp_assert(nkinds .LE. SIZE(rcut), cp_failure_level, &
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
       CALL fb_atomic_halo_nullify(halos(ihalo))
       CALL fb_atomic_halo_create(halos(ihalo), error)
    END DO
    CALL fb_atomic_halo_list_set(atomic_halos=atomic_halos, &
                                 nhalos=natoms_local, &
                                 max_nhalos=max_natoms_local, &
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
          CALL fb_atomic_halo_set(atomic_halo=halos(ihalo), &
                                  owner_id_in_halo=owner_id_in_halo, &
                                  natoms=nhalo_atoms, &
                                  halo_atoms=halo_atoms, &
                                  nelectrons=nelectrons_in_halo, &
                                  error=error)
          ! sort atomic halo
          CALL fb_atomic_halo_sort(halos(ihalo), error)
          CALL fb_atomic_halo_set(atomic_halo=halos(ihalo), &
                                  sorted=.TRUE., &
                                  error=error)
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



  SUBROUTINE fb_env_build_trial_fns_auto(fb_env, qs_env, nelectrons, maxocc, error)

    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    REAL(KIND=dp), INTENT(in) :: nelectrons
    REAL(KIND=dp), INTENT(in) :: maxocc
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_trial_fns_auto', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    TYPE(fb_trial_fns_obj) :: trial_fns
    TYPE(qs_kind_type), DIMENSION(:), POINTER :: qs_kind_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: basis_set_list
    TYPE(gto_basis_set_type), POINTER :: basis_set
    INTEGER :: basis_set_id, old_lshell, iset, nset, ikind, nkinds, &
               ishell, lshell, icgf, itrial, max_n_trial, counter, &
               ico, zz, stat
    INTEGER, DIMENSION(:), POINTER :: nfunctions, lmax, nshell
    INTEGER, DIMENSION(:,:), POINTER :: functions


    REAL(KIND=dp) :: zeff

!BEG:DEBUG:LT:2015/01/29
    INTEGER :: LT_UNIT
!END:DEBUG:LT:2015/01/29

!BEG:DEBUG:LT:2015/01/30
    LT_UNIT = 200
!END:DEBUG:LT:2015/01/30

!BEG:DEBUG:LT:2015/01/20
    write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: constructing trial functions"
!END:DEBUG:LT:2015/01/20

    failure = .FALSE.
    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)
    NULLIFY(nfunctions, functions, basis_set, basis_set_list, &
            atomic_kind_set, qs_kind_set)
    CALL fb_trial_fns_nullify(trial_fns)

    ! create a new trial_fn object
    CALL fb_trial_fns_create(trial_fns, error)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set, &
                    atomic_kind_set=atomic_kind_set, &
                    error=error)

    ! reading from the other parts of the code, it seemed that
    ! aux_fit_basis_set is only used when do_admm is TRUE. This can be
    ! seen from the calls to generate_qs_task_list subroutine in
    ! qs_create_task_list, found in qs_environment_methods.F:
    ! basis_set_id is only set as input parameter for do_admm
    ! calculations, and if not set, the task list is generated using
    ! the default basis_set_id=use_orb_basis_set.
    IF (qs_env%dft_control%do_admm) THEN
       basis_set_id = use_aux_fit_basis_set
    ELSE
       basis_set_id = use_orb_basis_set
    END IF

    nkinds = SIZE(qs_kind_set)

    ALLOCATE(nfunctions(nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    nfunctions = 0

    ALLOCATE(basis_set_list(nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL basis_set_list_setup(basis_set_list, basis_set_id, qs_kind_set)

    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)
       ! CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
       !                      z=zz)

       CALL get_qs_kind(qs_kind=qs_kind_set(ikind), &
                        zeff=zeff)

       bset1: DO iset = 1, nset
!          old_lshell = lmax(iset)
          old_lshell = -1
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                counter = counter + 1
                ! only include the first zeta orbitals
                IF ((lshell .GT. old_lshell) .AND. (counter .LE. nco(lshell))) THEN
                   nfunctions(ikind) = nfunctions(ikind) + 1
                END IF
             END DO
             ! we have got enough trial functions when we have enough
             ! basis functions to accomodate the number of electrons,
             ! AND that that we have included all the first zeta
             ! orbitals of an angular momentum quantum number l
             IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax(iset))) .AND. &
                 (maxocc * REAL(nfunctions(ikind),dp) .GE. zeff)) THEN
                EXIT bset1
             END IF
             old_lshell = lshell
          END DO
       END DO bset1
    END DO ! ikind

    ! now that we have the number of trial functions get the trial
    ! functions
    max_n_trial = MAXVAL(nfunctions)

!BEG:DEBUG:LT:2015/01/26
    write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: max_n_trial = ", max_n_trial
    write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: nkinds = ", nkinds
    write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: nfunctions = ", nfunctions
!END:DEBUG:LT:2015/01/26

    ALLOCATE(functions(max_n_trial,nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    functions(:,:) = 0
    ! redo the loops to get the trial function indices within the basis set
    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)

       ! CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
       !                      z=zz)

       CALL get_qs_kind(qs_kind=qs_kind_set(ikind), &
                        zeff=zeff)

!BEG:DEBUG:LT:2015/01/26
       write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: ikind = ", ikind
       write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: zeff = ", zeff
!END:DEBUG:LT:2015/01/26

       icgf = 0
       itrial = 0
       bset2: DO iset = 1, nset
!          old_lshell = lmax(iset)
          old_lshell = -1
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                icgf = icgf + 1
                counter = counter + 1
                ! only include the first zeta orbitals

!BEG:DEBUG:LT:2015/01/20
                write (LT_UNIT,*) &
                      "|LT| fb_env_build_trial_fns_auto: icg, counter, lshell, nco(lshell) = ", &
                      icgf, counter, lshell, nco(lshell)
!END:DEBUG:LT:2015/01/20

                IF ((lshell .GT. old_lshell) .AND. (counter .LE. nco(lshell))) THEN
                   itrial = itrial + 1
                   functions(itrial,ikind) = icgf
                END IF
             END DO
             ! we have got enough trial functions when we have more
             ! basis functions than the number of electrons (obtained
             ! from atomic z), AND that that we have included all the
             ! first zeta orbitals of an angular momentum quantum
             ! number l
             IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax(iset))) .AND. &
                 (maxocc * REAL(itrial,dp) .GE. zeff)) THEN
                EXIT bset2
             END IF
             old_lshell = lshell
          END DO
       END DO bset2
    END DO ! ikind

!BEG:DEBUG:LT:2015/01/20
    write (LT_UNIT,*) "|LT| fb_env_build_trial_fns_auto: trial_functions = ", functions
!END:DEBUG:LT:2015/01/20

    ! set trial_functions
    CALL fb_trial_fns_set(trial_fns=trial_fns, &
                          nfunctions=nfunctions, &
                          functions=functions, &
                          error=error)
    ! set fb_env
    CALL fb_env_set(fb_env=fb_env, &
                    trial_fns=trial_fns, &
                    error=error)

    ! cleanup
    DEALLOCATE(basis_set_list, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE fb_env_build_trial_fns_auto


  SUBROUTINE fb_dbcsr_copy_sparse_struct(matrix_out, matrix_in, error)

    TYPE(cp_dbcsr_type), INTENT(in) :: matrix_in
    TYPE(cp_dbcsr_type), INTENT(inout) :: matrix_out
    TYPE(cp_error_type), INTENT(INOUT) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_dbcsr_copy_sparse_struct', &
                                   routineP = moduleN//':'//routineN

    LOGICAL :: failure
    CHARACTER(LEN=1)                         :: symmetry
    TYPE(cp_dbcsr_iterator) :: iter
    INTEGER :: iatom, jatom, iblk, nblks, stat
    REAL(dp), DIMENSION(:, :), POINTER :: mat_block
    INTEGER, DIMENSION(:), ALLOCATABLE :: rows, cols

!BEG:DEBUG:LT:2015/02/02
    INTEGER :: LT_UNIT
!END:DEBUG:LT:2015/02/02

!BEG:DEBUG:LT:2015/02/02
    LT_UNIT = 200
!END:DEBUG:LT:2015/02/02

    failure = .FALSE.

    nblks = matrix_in%matrix%m%nblks
    ALLOCATE(rows(nblks), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(cols(nblks), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

!BEG:DEBUG:LT:2015/02/02
    write (LT_UNIT,*) "|LT| fb_dbcsr_copy_sparse_struct: nblks = ", nblks
    CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02

    iblk = 0
    CALL cp_dbcsr_iterator_start(iter, matrix_in)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, &
                                         iatom, &
                                         jatom, &
                                         mat_block, &
                                         iblk)

!BEG:DEBUG:LT:2015/02/02
       write (LT_UNIT,*) "|LT| fb_dbcsr_copy_sparse_struct: iblk = ", iblk
       write (LT_UNIT,*) "|LT| fb_dbcsr_copy_sparse_struct: rows(iblk) = ", iatom
       write (LT_UNIT,*) "|LT| fb_dbcsr_copy_sparse_struct: cols(iblk) = ", jatom
       CALL FLUSH(LT_UNIT)
!END:DEBUG:LT:2015/02/02

       rows(iblk) = iatom
       cols(iblk) = jatom
    END DO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_reserve_blocks(matrix_out, rows, cols, error=error)
    CALL cp_dbcsr_finalize(matrix_out, error=error)

    ! cleanup
    DEALLOCATE(rows, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(cols, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE fb_dbcsr_copy_sparse_struct





END MODULE qs_fb_env_methods
