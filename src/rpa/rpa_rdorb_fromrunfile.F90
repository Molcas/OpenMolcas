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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_RdOrb_FromRunfile()

use RPA_globals, only: CMO, EMO, l_CMO, l_EMO, l_OccEn, l_VirEn, nBas, nOcc, nOrb, nSym, nVir, OccEn, VirEn
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iUHF, iSym, i, nB, nB2, ip, ipO, ipV
logical(kind=iwp) :: Found
character(len=*), parameter :: SecNam = 'RPA_RdOrb_FromRunfile'
integer(kind=iwp), external :: RPA_iUHF

! Restricted (1) or unrestricted (2)
iUHF = RPA_iUHF()

! Allocate memory for CMO
l_CMO = nBas(1)*nOrb(1)
nB2 = nBas(1)**2
do iSym=2,nSym
  l_CMO = l_CMO+nBas(iSym)*nOrb(iSym)
  nB2 = nB2+nBas(iSym)**2
end do
call mma_allocate(CMO,l_CMO,iUHF,label='CMO(RPA)')

! Read CMO array(s) from Runfile
call Get_dArray_chk('Last orbitals',CMO(:,1),nB2)
if (iUHF == 2) then
  call Get_dArray('CMO_ab',CMO(:,2),nB2)
end if

! Allocate memory for orbital energies
nB = nBas(1)
do iSym=2,nSym
  nB = nB+nBas(iSym)
end do
do i=1,iUHF
  l_OccEn(i) = nOcc(1,i)
  l_VirEn(i) = nVir(1,i)
  do iSym=2,nSym
    l_OccEn(i) = l_OccEn(i)+nOcc(iSym,i)
    l_VirEn(i) = l_VirEn(i)+nVir(iSym,i)
  end do
end do
call mma_allocate(OccEn,maxval(l_OccEn),iUHF,label='OccEn')
call mma_allocate(VirEn,maxval(l_VirEn),iUHF,label='VirEn')

! Read orbital energies from Runfile
call qpg_dArray('OrbE',Found,l_EMO)
if ((.not. Found) .or. (l_EMO == 0)) then
  call RPA_Warn(3,SecNam//': Did not find OrbE')
end if
call mma_allocate(EMO,l_EMO,iUHF,label='OrbE')
call Get_dArray('OrbE',EMO(:,1),l_EMO)
if (l_EMO /= nB) then
  call RPA_Warn(3,SecNam//': unexpected EMO dimension')
end if
ip = 1
ipO = 1
ipV = 1
do iSym=1,nSym
  call dCopy_(nOcc(iSym,1),EMO(ip,1),1,OccEn(ipO,1),1)
  call dCopy_(nVir(iSym,1),EMO(ip+nOcc(iSym,1),1),1,VirEn(ipV,1),1)
  ip = ip+nOrb(iSym)
  ipO = ipO+nOcc(iSym,1)
  ipV = ipV+nVir(iSym,1)
end do
if (iUHF == 2) then
  call Get_dArray('OrbE_ab',EMO(:,2),l_EMO)
  ip = 1
  ipO = 1
  ipV = 1
  do iSym=1,nSym
    call dCopy_(nOcc(iSym,2),EMO(ip,2),1,OccEn(ipO,2),1)
    call dCopy_(nVir(iSym,2),EMO(ip+nOcc(iSym,2),2),1,VirEn(ipV,2),1)
    ip = ip+nOrb(iSym)
    ipO = ipO+nOcc(iSym,2)
    ipV = ipV+nVir(iSym,2)
  end do
end if

end subroutine RPA_RdOrb_FromRunfile
