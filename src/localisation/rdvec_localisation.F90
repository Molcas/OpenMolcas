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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!               2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine RdVec_Localisation(nSym,nBas,nOrb,IndT,CMO,Occ,EOrb,FName)
! Thomas Bondo Pedersen, July 2010.
!
! Read orbital info and return in a format suitable for module
! localisation: deleted orbitals are included (as zeros). This is a
! work-around to fix bugs when orbitals are deleted.

! nSym: number of irreps
! nBas: number of basis functions
! nOrb: number of orbitals
! IndT: type indices, dim: nBas
! CMO: MO coefficients, dim: nBas*nBas
! Occ: Occupation numbers, dim: nBas
! EOrb: dim: nBas
! FName: filename with input orbitals

#include "intent.fh"

use Localisation_globals, only: fileorb_id, isHDF5
#ifdef _HDF5_
use mh5, only: mh5_close_file
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp), intent(_OUT_) :: IndT(*)
real(kind=wp), intent(_OUT_) :: CMO(*), Occ(*), EOrb(*)
character(len=*), intent(in) :: FName
integer(kind=iwp) :: iErr, iSym, iUHF, iWarn, iWFType, k1, k2, l_CMO, Lu, n1, n2, nBasT, nOrbT
real(kind=wp) :: Dummy(1)
character(len=80) :: VTitle
integer(kind=iwp), allocatable :: Ind_(:)
real(kind=wp), allocatable :: CMO_(:), EOr_(:), Occ_(:)
character(len=*), parameter :: SecNam = 'RdVec_Localisation'

#include "warnings.h"

nBasT = nBas(1)
nOrbT = nOrb(1)
do iSym=2,nSym
  nBasT = nBasT+nBas(iSym)
  nOrbT = nOrbT+nOrb(iSym)
end do

l_CMO = nBas(1)*nOrb(1)
do iSym=2,nSym
  l_CMO = l_CMO+nBas(iSym)*nOrb(iSym)
end do
call mma_allocate(CMO_,l_CMO,label='CMO_')
call mma_allocate(Occ_,nOrbT,label='Occ_')
call mma_allocate(EOr_,nOrbT,label='EOr_')
call mma_allocate(Ind_,nBasT,label='Ind_')

if (isHDF5) Then
  call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO_,Occ_,EOr_,Ind_)
# ifdef _HDF5_
  call mh5_close_file(fileorb_id)
# endif
else
  Lu = 75
  iUHF = 0  ! restricted HF
  iWarn = 2 ! abend if nBas/nOrb info is inconsistent
  iErr = -1 ! init return code
  iWFType = -1 ! init wave function type
  Dummy(1) = huge(Dummy) ! dummy variable
  call RdVec_(FName,Lu,'COEI',iUHF,nSym,nBas,nOrb,CMO_,Dummy,Occ_,Dummy,EOr_,Dummy,Ind_,VTitle,iWarn,iErr,iWFType)
  if (iErr /= 0) then
    call WarningMessage(2,SecNam//': Non-zero return code from RdVec_')
    write(u6,'(A,A,I9)') SecNam,': RdVec_ returned code',iErr
    call xFlush(u6)
    call xQuit(_RC_IO_ERROR_READ_)
  end if
  write(u6,*)
  write(u6,'(A)') ' Header from vector file:'
  write(u6,*)
  write(u6,'(A)') trim(VTitle)
end if
write(u6,*)

k1 = 0
k2 = 0
do iSym=1,nSym
  n1 = nBas(iSym)*nOrb(iSym)
  n2 = nBas(iSym)*nBas(iSym)
  CMO(k2+1:k2+n1) = CMO_(k1+1:k1+n1)
  CMO(k2+n1+1:k2+n2) = Zero
  k1 = k1+n1
  k2 = k2+n2
end do

k1 = 0
k2 = 0
do iSym=1,nSym
  n1 = nOrb(iSym)
  n2 = nBas(iSym)
  Occ(k2+1:k2+n1) = Occ_(k1+1:k1+n1)
  Occ(k2+n1+1:k2+n2) = Zero
  EOrb(k2+1:k2+n1) = EOr_(k1+1:k1+n1)
  EOrb(k2+n1+1:k2+n2) = huge(Dummy)
  IndT(k2+1:k2+n1) = Ind_(k1+1:k1+n1)
  IndT(k2+n1+1:k2+n2) = 7
  k1 = k1+n1
  k2 = k2+n1
end do

call mma_deallocate(CMO_)
call mma_deallocate(Occ_)
call mma_deallocate(EOr_)
call mma_deallocate(Ind_)

end subroutine RdVec_Localisation
