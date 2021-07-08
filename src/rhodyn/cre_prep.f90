subroutine cre_prep

  use rhodyn_data
  use mh5
  implicit none

! creating intermediate preparation file
  prep_id = mh5_create_file('RDPREP')

  call mh5_init_attr (prep_id,'MOLCAS_MODULE','RHODYN')

! CI coefficients
  prep_ci = mh5_create_dset_real (prep_id, &
       'CI_COEFF', 3, [maxnconf,maxlroots,N])
  call mh5_init_attr(prep_ci, 'description','CI coefficients')

! SF hamiltonians
  prep_hcsf = mh5_create_dset_real (prep_id, &
       'SFS_HAM', 3, [maxnconf,maxnconf,N])
  call mh5_init_attr(prep_hcsf, 'description','SF Hamiltonians')

! U_CI
  prep_uci = mh5_create_dset_real (prep_id, &
       'U_CI', 2, [nconftot,lrootstot])
  call mh5_init_attr(prep_uci, 'description', &
       'trafo matrix accounting for spin-degeneracy')

! SO-Hamiltonian in CSF basis, V_CSF
  prep_vcsfr = mh5_create_dset_real (prep_id, &
        'V_CSF_R', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_vcsfr, 'description', &
        'SO-Hamiltonian in CSF basis, real part')
  prep_vcsfi = mh5_create_dset_real (prep_id, &
        'V_CSF_I', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_vcsfi, 'description', &
        'SO-Hamiltonian in CSF basis, imaginary part')

! Full Hamiltonian in CSF basis, HTOT_CSF
  prep_fhr = mh5_create_dset_real (prep_id, &
        'FULL_H_R', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_fhr, 'description', &
        'SO-Hamiltonian in CSF basis, real part')
  prep_fhi = mh5_create_dset_real (prep_id, &
        'FULL_H_I', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_fhi, 'description', &
        'SO-Hamiltonian in CSF basis, imaginary part')

! CSF2SO
  prep_csfsor = mh5_create_dset_real (prep_id, &
        'CSF2SO_R', 2, [nconftot,lrootstot])
  call mh5_init_attr(prep_csfsor, 'description','CSF2SO_R')
  prep_csfsoi = mh5_create_dset_real (prep_id, &
        'CSF2SO_I', 2, [nconftot,lrootstot])
  call mh5_init_attr(prep_csfsoi, 'description','CSF2SO_I')

! Dipole SO matrices
  prep_dipoler = mh5_create_dset_real (prep_id, &
        'DIP_MOM_R', 3, [lrootstot,lrootstot,3])
  call mh5_init_attr(prep_dipoler, 'description', &
        'Dipole SO matrices, real part')
  prep_dipolei = mh5_create_dset_real (prep_id, &
        'DIP_MOM_I', 3, [lrootstot,lrootstot,3])
  call mh5_init_attr(prep_dipolei, 'description', &
        'Dipole SO matrices, imaginary part')

! Initial density matrix
  prep_dm_r = mh5_create_dset_real (prep_id, &
        'DM0_R', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_dm_r, 'description', &
        'Initial density matrix, real part')
  prep_dm_i = mh5_create_dset_real (prep_id, &
        'DM0_I', 2, [nconftot,nconftot])
  call mh5_init_attr(prep_dm_i, 'description', &
        'Initial density matrix, imaginary part')

! Dyson amplitudes matrix
  prep_do = mh5_create_dset_real(prep_id,'DYSAMP', &
         2,[lrootstot,lrootstot])
  call mh5_init_attr(prep_do, 'description', &
        'Matrix of Dyson amplitudes in SF basis if'// &
        'number of spin manifolds >2, otw in SO')

end
