!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************

subroutine cre_prep()

use rhodyn_data, only: flag_dyson, flag_so, lrootstot, maxlroots, maxnconf, N, nconftot, prep_ci, prep_csfsoi, prep_csfsor, &
                       prep_dipolei, prep_dipoler, prep_dm_i, prep_dm_r, prep_do, prep_fhi, prep_fhr, prep_hcsf, prep_id, &
                       prep_uci, prep_vcsfi, prep_vcsfr
use mh5, only: mh5_create_dset_real, mh5_create_file, mh5_init_attr

implicit none

! creating intermediate preparation file
prep_id = mh5_create_file('RDPREP')

call mh5_init_attr(prep_id,'MOLCAS_MODULE','RHODYN')

! CI coefficients
prep_ci = mh5_create_dset_real(prep_id,'CI_COEFF',3,[maxnconf,maxlroots,N])
call mh5_init_attr(prep_ci,'description','CI coefficients')

! SF hamiltonians
prep_hcsf = mh5_create_dset_real(prep_id,'SFS_HAM',3,[maxnconf,maxnconf,N])
call mh5_init_attr(prep_hcsf,'description','SF Hamiltonians')

! U_CI
prep_uci = mh5_create_dset_real(prep_id,'U_CI',2,[nconftot,lrootstot])
call mh5_init_attr(prep_uci,'description','trafo matrix accounting for spin-degeneracy')

! SO-Hamiltonian in CSF basis, V_CSF
if (flag_so) then
  prep_vcsfr = mh5_create_dset_real(prep_id,'V_CSF_R',2,[nconftot,nconftot])
  call mh5_init_attr(prep_vcsfr,'description','SO-Hamiltonian in CSF basis, real part')
  prep_vcsfi = mh5_create_dset_real(prep_id,'V_CSF_I',2,[nconftot,nconftot])
  call mh5_init_attr(prep_vcsfi,'description','SO-Hamiltonian in CSF basis, imaginary part')
end if

! Full Hamiltonian in CSF basis, HTOT_CSF
prep_fhr = mh5_create_dset_real(prep_id,'FULL_H_R',2,[nconftot,nconftot])
call mh5_init_attr(prep_fhr,'description','Hamiltonian in CSF basis, real part')
prep_fhi = mh5_create_dset_real(prep_id,'FULL_H_I',2,[nconftot,nconftot])
call mh5_init_attr(prep_fhi,'description','Hamiltonian in CSF basis, imaginary part')

! CSF2SO
if (flag_so) then
  prep_csfsor = mh5_create_dset_real(prep_id,'CSF2SO_R',2,[nconftot,lrootstot])
  call mh5_init_attr(prep_csfsor,'description','CSF2SO_R')
  prep_csfsoi = mh5_create_dset_real(prep_id,'CSF2SO_I',2,[nconftot,lrootstot])
  call mh5_init_attr(prep_csfsoi,'description','CSF2SO_I')
end if

! Dipole matrices
prep_dipoler = mh5_create_dset_real(prep_id,'DIP_MOM_R',3,[lrootstot,lrootstot,3])
call mh5_init_attr(prep_dipoler,'description','Dipole matrices, real part')
prep_dipolei = mh5_create_dset_real(prep_id,'DIP_MOM_I',3,[lrootstot,lrootstot,3])
call mh5_init_attr(prep_dipolei,'description','Dipole matrices, imaginary part')

! Initial density matrix in CSF basis
prep_dm_r = mh5_create_dset_real(prep_id,'DM0_R',2,[nconftot,nconftot])
call mh5_init_attr(prep_dm_r,'description','Initial density matrix, real part')
prep_dm_i = mh5_create_dset_real(prep_id,'DM0_I',2,[nconftot,nconftot])
call mh5_init_attr(prep_dm_i,'description','Initial density matrix, imaginary part')

! Dyson amplitudes matrix
if (flag_dyson) then
  prep_do = mh5_create_dset_real(prep_id,'DYSAMP',2,[lrootstot,lrootstot])
  call mh5_init_attr(prep_do,'description','Matrix of Dyson amplitudes in SF basis if number of spin manifolds >2, otw in SO')
end if

end subroutine cre_prep
