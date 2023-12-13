!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine cre_locwfn()

! Create a wavefunction file.
! If another .local.h5 file already exists, it will be overwritten.
#ifdef _HDF5_
use Localisation_globals, only: nBas, nSym, wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
use mh5, only: mh5_create_file, mh5_init_attr, mh5_create_dset_real, mh5_create_dset_str
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nBasTot, nSqrTot

! create a new wavefunction file!
wfn_fileid = mh5_create_file('LOCWFN')

! set module type
call mh5_init_attr(wfn_fileid,'MOLCAS_MODULE','LOCALISATION')

! copy basic molecular information to the HDF5 file
call run2h5_molinfo(wfn_fileid)
call one2h5_ovlmat(wfn_fileid,nsym,nbas)
call one2h5_fckint(wfn_fileid,nsym,nbas)

nBasTot = sum(nBas(1:nSym))
nSqrTot = sum(nBas(1:nSym)**2)

! typestring
wfn_tpidx = mh5_create_dset_str(wfn_fileid,'MO_TYPEINDICES',1,[nBasTot],1)
call mh5_init_attr(wfn_tpidx,'DESCRIPTION', &
                   'Type index of the molecular orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')
! molecular orbital coefficients
wfn_mocoef = mh5_create_dset_real(wfn_fileid,'MO_VECTORS',1,[nSqrTot])
call mh5_init_attr(wfn_mocoef,'DESCRIPTION', &
                   'Coefficients of the molecular orbitals, arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
! molecular orbital occupation numbers
wfn_occnum = mh5_create_dset_real(wfn_fileid,'MO_OCCUPATIONS',1,[nBasTot])
call mh5_init_attr(wfn_occnum,'DESCRIPTION', &
                   'Occupation numbers of the molecular orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')
! molecular orbital energies
wfn_orbene = mh5_create_dset_real(wfn_fileid,'MO_ENERGIES',1,[nBasTot])
call mh5_init_attr(wfn_orbene,'DESCRIPTION', &
                   'Orbital energies of the molecular orbitals arranged as blocks of size [NBAS(i)], i=1,#irreps')
#endif

end subroutine cre_locwfn

!-----------------------------------------------------------------------

subroutine cls_locwfn()

#ifdef _HDF5_
use Localisation_globals, only: wfn_fileid
use mh5, only: mh5_close_file

implicit none

call mh5_close_file(wfn_fileid)
#endif

end subroutine cls_locwfn
