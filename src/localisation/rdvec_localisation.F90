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
!***********************************************************************

subroutine RdVec_Localisation(nSym,nBas,nOrb,IndT,CMO,Occ,EOrb,FName)

! Thomas Bondo Pedersen, July 2010.
!
! Read orbital info and return in a format suitable for module
! localisation: deleted orbitals are included (as zeros). This is a
! work-around to fix bugs when orbitals are deleted.

implicit none
integer nSym ! number of irreps
integer nBas(nSym)  ! number of basis functions
integer nOrb(nSym)  ! number of orbitals
integer IndT(*)  ! type indices, dim: nBas
real*8 CMO(*)   ! MO coefficients, dim: nBas*nBas
real*8 Occ(*)   ! Occupation numbers, dim: nBas
real*8 EOrb(*)  ! EOrb, dim: nBas
character*(*) FName ! filename with input orbitals
#include "warnings.fh"
#include "WrkSpc.fh"

character*18 SecNam
parameter(SecNam='RdVec_Localisation')

character*80 VTitle

integer nBasT, nOrbT, iSym
integer ip_CMO, l_CMO
integer ip_Occ, l_Occ
integer ip_EOr, l_EOr
integer ip_Ind, l_Ind
integer Lu, iUHF, iWarn, iErr, iWFType, k1, k2, i
integer mylen
external mylen

real*8 Dummy(1)

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
l_Occ = nOrbT
l_EOr = nOrbT
l_Ind = nBasT
call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)

Lu = 75
iUHF = 0  ! restricted HF
iWarn = 2 ! abend if nBas/nOrb info is inconsistent
iErr = -1 ! init return code
iWFType = -1 ! init wave function type
Dummy(1) = 9.9d9 ! dummy variable
call RdVec_(FName,Lu,'COEI',iUHF,nSym,nBas,nOrb,Work(ip_CMO),Dummy,Work(ip_Occ),Dummy,Work(ip_EOr),Dummy,iWork(ip_Ind),VTitle, &
            iWarn,iErr,iWFType)
if (iErr /= 0) then
  call WarningMessage(2,SecNam//': Non-zero return code from RdVec_')
  write(6,'(A,A,I9)') SecNam,': RdVec_ returned code',iErr
  call xFlush(6)
  call xQuit(_RC_IO_ERROR_READ_)
end if
write(6,*)
write(6,'(A)') ' Header from vector file:'
write(6,*)
write(6,'(A)') VTitle(:mylen(VTitle))
write(6,*)

k1 = ip_CMO
k2 = 1
do iSym=1,nSym
  call dCopy_(nBas(iSym)*nOrb(iSym),Work(k1),1,CMO(k2),1)
  call Cho_dZero(CMO(k2+nBas(iSym)*nOrb(iSym)),nBas(iSym)*(nBas(iSym)-nOrb(iSym)))
  k1 = k1+nBas(iSym)*nOrb(iSym)
  k2 = k2+nBas(iSym)*nBas(iSym)
end do

k1 = ip_Occ
k2 = 1
do iSym=1,nSym
  call dCopy_(nOrb(iSym),Work(k1),1,Occ(k2),1)
  call Cho_dZero(Occ(k2+nOrb(iSym)),nBas(iSym)-nOrb(iSym))
  k1 = k1+nOrb(iSym)
  k2 = k2+nBas(iSym)
end do

k1 = ip_EOr
k2 = 1
Dummy(1) = 9.9d9
do iSym=1,nSym
  call dCopy_(nOrb(iSym),Work(k1),1,EOrb(k2),1)
  call dCopy_(nBas(iSym)-nOrb(iSym),Dummy(1),0,EOrb(k2+nOrb(iSym)),1)
  k1 = k1+nOrb(iSym)
  k2 = k2+nBas(iSym)
end do

k1 = ip_Ind
k2 = 1
do iSym=1,nSym
  call iCopy(nOrb(iSym),iWork(k1),1,IndT(k2),1)
  do i=nOrb(iSym),nBas(iSym)-1
    IndT(k2+i) = 7
  end do
  k1 = k1+nOrb(iSym)
  k2 = k2+nBas(iSym)
end do

call GetMem('Ind_','Free','Inte',ip_Ind,l_Ind)
call GetMem('EOr_','Free','Real',ip_EOr,l_EOr)
call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

end subroutine RdVec_Localisation
