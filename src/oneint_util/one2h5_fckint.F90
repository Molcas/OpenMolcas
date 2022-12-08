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

#include "compiler_features.h"
#ifdef _HDF5_

subroutine one2h5_fckint(fileid,nSym,nBas)
! IFG: read atomic Fock matrix from the 1-electron integral file
! and write it to the HDF5 file specified with fileid.
! This routine does nothing if HDF5 is not supported.
!
! Datasets:
!   AO_FOCKINT_MATRIX

use OneDat, only: sNoNuc, sNoOri
use mh5, only: mh5_close_dset, mh5_create_dset_real, mh5_init_attr, mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: fileid, nSym, nBas(nSym)
integer(kind=iwp) :: dsetid, iComp, iOff1, iOff2, iOpt, iRc, iSyLbl, iSym, nb, nbast, nbast1, nbast2
character(len=8) :: Label
real(kind=wp), allocatable :: FAO(:), Scr(:,:)

nbast = 0
nbast1 = 0
nbast2 = 0
do iSym=1,nSym
  nb = nBas(iSym)
  nbast = nbast+nb
  nbast1 = nbast1+(nb*(nb+1))/2
  nbast2 = nbast2+nb**2
end do

! atomic orbital Fock matrix
dsetid = mh5_create_dset_real(fileid,'AO_FOCKINT_MATRIX',1,[NBAST2])
call mh5_init_attr(dsetid,'DESCRIPTION','Fock matrix of the atomic orbitals, arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

call mma_allocate(FAO,NBAST1)
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'FckInt  '
call RdOne(iRc,iOpt,Label,iComp,FAO,iSyLbl)
iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  nb = nBas(iSym)
  if (nb > 0) then
    call mma_allocate(Scr,nb,nb,label='Scr')
    call Square(FAO(1+iOff1),Scr,1,nb,nb)
    call mh5_put_dset(dsetid,Scr,[nb*nb],[iOff2])
    call mma_deallocate(Scr)
  end if
  iOff1 = iOff1+nb*(nb+1)/2
  iOff2 = iOff2+nb**2
end do
call mma_deallocate(FAO)

call mh5_close_dset(dsetid)

end subroutine one2h5_fckint

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(one2h5_fckint)

#endif
