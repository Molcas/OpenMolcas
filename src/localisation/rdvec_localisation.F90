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

! nSym: number of irreps
! nBas: number of basis functions
! nOrb: number of orbitals
! IndT: type indices, dim: nBas
! CMO: MO coefficients, dim: nBas*nBas
! Occ: Occupation numbers, dim: nBas
! EOrb: dim: nBas
! FName: filename with input orbitals

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp), intent(out) :: IndT(*)
real(kind=wp), intent(out) :: CMO(*), Occ(*), EOrb(*)
character(len=*), intent(in) :: FName
#include "warnings.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iErr, ip_CMO, ip_EOr, ip_Ind, ip_Occ, iSym, iUHF, iWarn, iWFType, k1, k2, l_CMO, l_EOr, l_Ind, l_Occ, Lu, &
                     nBasT, nOrbT
real(kind=wp) :: Dummy(1)
character(len=80) :: VTitle
character(len=18), parameter :: SecNam = 'RdVec_Localisation'

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
Dummy(1) = huge(Dummy) ! dummy variable
call RdVec_(FName,Lu,'COEI',iUHF,nSym,nBas,nOrb,Work(ip_CMO),Dummy,Work(ip_Occ),Dummy,Work(ip_EOr),Dummy,iWork(ip_Ind),VTitle, &
            iWarn,iErr,iWFType)
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
write(u6,*)

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
Dummy(1) = huge(Dummy)
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
