MODULE qs_fb_env_methods

  USE kinds, ONLY: dp, &
                   default_string_length
  USE qs_fb_env_types, ONLY: fb_env_obj, &
                             fb_env_get, &
                             fb_env_set
  USE qs_environment_types,    ONLY: qs_environment_type, &
                                     get_qs_env
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
                                     fb_atomic_halo_list_has_data, &
                                     fb_atomic_halo_list_nullify
  USE particle_types,        ONLY: particle_type
  USE atomic_kind_types,     ONLY: atomic_kind_type, &
                                   get_atomic_kind
  USE qs_neighbor_list,      ONLY: pair_radius_setup
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
                                cp_dbcsr_release
  USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                 cp_dbcsr_sm_fm_multiply
  USE input_constants, ONLY: use_aux_fit_basis_set,&
                             use_orb_basis_set, &
                             cholesky_inverse
  USE orbital_pointers, ONLY: ncoset, &
                              nco
  USE qs_mo_types, ONLY: mo_set_p_type, &
                         get_mo_set, &
                         init_mo_set, &
                         set_mo_occupation, &
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
                         cp_fm_release
  USE qs_fb_trial_fns_types, ONLY: fb_trial_fns_obj, &
                                   fb_trial_fns_nullify, &
                                   fb_trial_fns_get
  USE array_types, ONLY: array_i1d_obj, &
                         array_nullify
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


  SUBROUTINE fb_env_do_diag(fb_env, qs_env, matrix_ks, matrix_s, &
                            diis_step, error)
    TYPE(fb_env_obj), INTENT(in) :: fb_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_s
    LOGICAL, INTENT(INOUT) :: diis_step
    TYPE(cp_error_type), INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_do_diag', &
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
    TYPE(cp_fm_type), OINTER :: mo_coeff_filtered, mo_coeff, &
                                fm_matrix_filtered_s, fm_matrix_filtered_ks, &
                                fm_matrix_work
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(section_vals_type), POINTER :: scf_section
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos, mos_filtered
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_filter
    TYPE(cp_dbcsr_type) :: matrix_tmp, matrix_filtered_ks, &
                           matrix_filtered_p, matrix_filtered_s
    TYPE(fb_atomic_halo_list_obj) :: atomic_halos
    TYPE(fb_trial_fns_obj) :: trial_fns
    CHARACTER(len=2) :: spin_string
    CHARACTER(len=default_string_length) :: name
    TYPE(array_i1d_obj) :: filtered_rowORcol_block_sizes, &
                           original_rowORcol_block_sizes
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(qs_matrix_pools_type), POINTER :: my_mpools
    TYPE(cp_fm_struct_type), POINTER  :: fm_struct


    failure = .FALSE.

    NULLIFY(scf_env, scf_control, scf_section, para_env, blacs_env, &
            mos, mos_filtered, matrix_filter, particle_set)
    CALL fb_atomic_halo_list_nullify(atomic_halos)
    CALL fb_trial_fns_nullify(trial_fns)
    CALL array_nullify(original_rowORcol_block_sizes)
    CALL array_nullify(filtered_rowORcol_block_sizes)

    nspin = SIZE(matrix_ks)

    ! get qs_env information
    CALL get_qs_env(qs_env=qs_env, &
                    scf_env=scf_env, &
                    scf_control=scf_control, &
                    scf_senction=scf_section, &
                    para_env=para_env, &
                    blacs_env=blacs_env
                    particle_set=particle_set, &
                    mos=mos)

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
       ! get filter matrix name
       WRITE (spin_string, FMT="(I1)"), ispin
       name = TRIM("FILTER MATRIX SPIN "//spin_string)
       CALL compress(name)
       CALL uppercase(name)
       ! calculate filter matrix (matrix_s(1) is the overlap, the rest
       ! in the array are its derivatives)
       CALL fb_fltrmat_build(H_mat=matrix_ks(ispin)%matrix, &
                             S_mat=matrix_s(1)%matrix, &
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

    ! get the filtered density matrix and then convert back to the
    ! full basis version in scf_env ready to be used outside this
    ! subroutine
    CALL cp_dbcsr_init(matrix_filtered_p, error=error)
    CALL cp_dbcsr_create(matrix=matrix_filtered_p, &
                         name=TRIM("FILTERED_MATRIX_P"), &
                         ! all DBCSR matrices in CP2K uses the same distribution
                         dist=cp_dbcsr_distribution(matrix_s(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_p, error=error)

    DO ispin = nspin
       ! calculate matrix_filtered_p
       CALL calculate_density_matrix(mos_filtered(ispin)%mo_set, &
                                     matrix_filtered_p, &
                                     error=error)
       ! convert back to full basis p
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_filtered_p, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                              matrix_tmp, matrix_filter(ispin)%matrix, &
                              0.0_dp, scf_env%p_mix_new(ispin)%matrix, error=error)
    END DO

    ! release temporary matrices
    CALL cp_dbcsr_release(matrix_tmp, error)
    CALL cp_dbcsr_release(matrix_filtered_p, error)

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

    NULLIFY(eigenvalues_filtered, occ_filtered, mo_coeff_filtered, &
            eigenvalues, occ, mo_coeff)
    DO ispin = nspin
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
                       mu=mu_filtered)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       occupation_numbers=occ, &
                       eigenvalues=eigenvalues, &
                       mo_coeff=mo_coeff)
       eigenvalues(:) = 0.0_dp
       eigenvalues(1:nmo_filtered) = eigenvalues_filtered(1:nmo_filtered)
       occ(:) = 0.0_dp
       occ(1:nmo_filtered) = occ_filtered(1:nmo_filtered)
       CALL cp_dbcsr_sm_fm_multiply(matrix_filter(ispin)%matrix, &
                                    mo_coeff_filtered, &
                                    mo_coeff, &
                                    nmo_filtered, &
                                    error=error)
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
    REAK(KIND=dp) :: kind_radius
    INTEGER :: ikind, basis_set_id, nkinds, stat

    NULLIFY(rcut, qs_kind_set)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set)

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

    ALLOCATE(basis_set_list(nkinds), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL basis_set_list_setup(basis_set_list, basis_set_id, qs_kind_set)

    DO ikind = 1, nkinds
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              kind_radius=kind_radius)
       rcut(ikind) = kind_radius
    END DO

    CALL fb_env_set(fb_env=fb_env, &
                    rcut=rcut, &
                    error=error)

    ! cleanup
    DEALLOCATE(basis_set_list, STAT=stat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)

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
                    natoms=natoms_global, &
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
          CALL fb_atomic_halo_set(halos(ihalo) &
                                  owner_id_in_halo=owner_id_in_halo, &
                                  natoms=nhalo_atoms, &
                                  halo_atoms=halo_atoms, &
                                  nelectrons=nelectrons_in_halo, &
                                  error=error)
          ! sort atomic halo
          CALL fb_atomic_halo_sort(halos(ihalo), error)

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
    TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: basis_set_list
    TYPE(gto_basis_set_type), POINTER :: basis_set
    INTEGER :: basis_set_id, old_lshell, iset, nset, ikind, nkinds, &
               ishell, lshell, icgf, itrial, max_n_trial, counter, &
               ico, stat
    INTEGER, DIMENSION(:), POINTER :: nfunctions, lmax, nshell
    INTEGER, DIMENSION(:,:), POINTER :: functions


    failure = .FALSE.
    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)
    NULLIFY(nfunctions, functions, basis_set, basis_set_list, &
            qs_kind_set, atomic_kind_set)
    CALL fb_trial_fns_nullify(trial_fns)

    ! create a new trial_fn object
    CALL fb_trial_fns_create(trial_fns, error)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set)

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

    ALLOCATE(basis_set_list(nkinds), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL basis_set_list_setup(basis_set_list, basis_set_id, qs_kind_set)

    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)
       bset: DO iset = 1, nset
          old_lshell = lmax(iset)
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                counter = counter + 1
                ! only include the first zeta orbitals
                IF (counter .LE. nco(lshell)) THEN
                   nfunctions(ikind) = nfunctions(ikind) + 1
                END IF
             END DO
             ! we have got enough trial functions when we have enough
             ! basis functions to accomodate the number of electrons,
             ! AND that that we have included all the first zeta
             ! orbitals of an angular momentum quantum number l
             IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax(iset))) .AND. &
                 (maxocc * REAL(nfunctions(ikind),dp) .GE. nelectrons)) THEN
                EXIT bset
             END IF
             old_lshell = lshell
          END DO
       END DO bset
    END DO ! ikind

    ! now that we have the number of trial functions get the trial
    ! functions
    max_n_trial = MAXVAL(nfunctions)
    ALLOCATE(functions(max_n_trial,nkinds), STAT=stat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    functions(:,:) = 0
    ! redo the loops to get the trial function indices within the basis set
    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=gto_basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)
       icgf = 0
       bset: DO iset = 1, nset
          old_lshell = lmax(iset)
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                icgf = icgf + 1
                counter = counter + 1
                ! only include the first zeta orbitals
                IF (counter .LE. nco(lshell)) THEN
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
                 (maxocc * REAL(nfunctions(ikind),dp) .GE. nelectrons)) THEN
                EXIT bset
             END IF
             old_lshell = lshell
          END DO
       END DO bset
    END DO ! ikind

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
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE fb_env_build_trial_fns_auto



  ! SUBROUTINE fb_env_build_filter_matrix(fb_env, qs_env, error)
  !   TYPE(fb_env_obj), POINTER :: fb_env
  !   TYPE(qs_environment_type), POINTER :: qs_env
  !   TYPE(cp_error_type), INTENT(INOUT) :: error

  !   CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_env_build_filter_matrix', &
  !                                  routineP = moduleN//':'//routineN

  !   LOGICAL :: failure

  !   TYPE(cp_para_env_type), POINTER :: para_env
  !   TYPE(dbcsr_distribution_obj), POINTER :: dbcsr_dist
  !   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mo_list
  !   TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
  !   TYPE(atomic_kind_type), POINTER :: atomic_kind
  !   TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: mat_S, mat_H
  !   TYPE(fb_atomic_halo_list_obj) :: atomic_halos
  !   TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
  !   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: atomic_H, atomic_S
  !   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: filtered_eigenvectors, &
  !                                                 atomic_filter_matrix
  !   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: filter_function, &
  !                                               eigenvalues
  !   REAL(KIND=dp) :: filter_temperature, fermi_level, maxocc
  !   INTEGER :: ispin, nblkrows_total, nblkcols_total, stat, jhalo, nhalos, &
  !              this_atom, this_atom_in_halo, ii, jj, mat_dim, nrows, ncols, &
  !              natoms_in_halo, ihalo_atom, jhalo_atom, iatom_global, jatom_global, &
  !              natoms_global, ikind, jkind, nspin
  !   INTEGER, DIMENSION(:), POINTER :: halo_atoms
  !   INTEGER, DIMENSION(:), ALLOCATABLE :: atomic_H_blk_row_start, &
  !                                         atomic_H_blk_col_start, &
  !                                         atomic_S_blk_row_start, &
  !                                         atomic_S_blk_col_start
  !   INTEGER, DIMENSION(:), POINTER :: cbs
  !   TYPE(array_i1d_obj) :: col_blk_size, row_blk_size

  !   failure = .FALSE.

  !   CPPostcondition(fb_env_has_data(fb_env), &
  !                   cp_failure_level, routineP, error, failure)

  !   NULLIFY(para_env, mo_list, particle_set, atomic_kind, mat_S, mat_H, halos, &
  !           atomic_H, atomic_S, atomic_H_blk_row_start, atomic_H_blk_col_start, &
  !           atomic_S_blk_row_start, atomic_S_blk_col_start, cbs)
  !   CALL array_nullify(row_blk_size)
  !   CALL array_nullify(col_blk_size)
  !   CALL fb_atomic_halo_list_nullify(atomic_halos)

  !   CALL fb_env_get(fb_env=fb_env, &
  !                   filter_temperature=filter_temperature, &
  !                   error=error)

  !   CALL get_qs_env(qs_env=qs_env, &
  !                   para_env=para_env, &
  !                   mos=mo_list, &
  !                   particle_set=particle_set, &
  !                   matrix_s=mat_S, &
  !                   matrix_h=mat_H, &
  !                   dbcsr_dist=dbcsr_dist, &
  !                   error=error)

  !   CPPrecondition(SIZE(mat_H)==SIZE(mat_S), cp_failure_level, routineP, error, failure)
  !   nspin = SIZE(mat_H)

  !   ! get trial function information
  !   CALL fb_trial_fns_get(trial_fns=trial_fns, &
  !                         nfunctions=ntfns, &
  !                         functions=tfns, &
  !                         error=error)

  !   ! allocate the row_blk_sizes and col_blk_sizes arrays for
  !   ! constructing the filter matrix in DBCSR format
  !   natoms_global = SIZE(particle_set)
  !   ALLOCATE(cbs(natoms_global), STAT=stat)
  !   CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !   rbs = 0
  !   cbs = 0
  !   ! row_blk_size is same as S or H matrix
  !   CALL cp_dbcsr_get_info(mat_S(ispin=1)%matrix, &
  !                          nblkrows_total=nblkrows_total, &
  !                          nblkcols_total=nblkcols_total, &
  !                          row_blk_size=row_blk_size)
  !   ! work out col_blk_size
  !   DO iatom_global = 1, natoms_global
  !      atomic_kind => particle_set(iatom_global)%atomic_kind
  !      CALL get_atomic_kind(atomic_kind=atomic_kind, &
  !                           kind_number=ikind)
  !      cbs(iatom_global) = ntfns(ikind)
  !   END DO
  !   ! do not deallocate cbs if gift=.TRUE. as col_blk_sizes will only point to cbs
  !   CALL array_new(col_blk_sizes, cbs, gift=.TRUE.)
  !   ! the filter matrix is non-square and is always non-symmetric
  !   symmetry_string = dbcsr_type_no_symmetry
  !   ! create empty filter matrices for each spin
  !   DO ispin = 1, nspin
  !      WRITE (spin_string, FMT="(I1)"), ispin
  !      name = TRIM("FILTER MATRIX SPIN"//spin_string)
  !      CALL compress(name)
  !      CALL uppercase(name)
  !      ALLOCATE(filter_matrix(ispin)%matrix, STAT=stat)
  !      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
  !      CALL cp_dbcsr_init(filter_matrix(ispin)%matrix, error=error)
  !      CALL cp_dbcsr_create(matrix=filter_matrix(ispin)%matrix, &
  !                           name=name, &
  !                           dist=dbcsr_dist, &
  !                           matrix_type=symmetry_string, &
  !                           row_blk_size=row_blk_size, &
  !                           col_blk_size=col_blk_size, &
  !                           nblks=0, &
  !                           nze=0, &
  !                           error=error)
  !   END DO
  !   ! generate atomic halos
  !   CALL fb_env_build_atomic_halos(fb_env, error)
  !   CALL fb_env_get(fb_env=fb_env, &
  !                   atomic_halos=atomic_halos, &
  !                   error=error)
  !   CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
  !                                nhalos=nhalos, &
  !                                halos=halos, &
  !                                error=error)
  !   DO ispin = 1, nspin
  !      DO jhalo = 1, nhalos
  !         ! construct the filter matrix block by block
  !         CALL fb_fltrmat_add_blkcol(mat_H(ispin)%matrix, &
  !                                    mat_S(ispin)%matrix, &
  !                                    halo(jhalo), &
  !                                    trial_fns, &
  !                                    para_env, &
  !                                    mo_list(ispin)%mo_set, &
  !                                    filter_temperature, &
  !                                    filter_matrix(ispin)%matrix, &
  !                                    error)
  !      END DO  ! ihalo
  !      ! finalise the filter matrix
  !      CALL dbcsr_finalize(filter_matrix(ispin)%matrix)
  !   END DO ! ispin

  END SUBROUTINE fb_env_build_filter_matrix




END MODULE qs_fb_env_methods
