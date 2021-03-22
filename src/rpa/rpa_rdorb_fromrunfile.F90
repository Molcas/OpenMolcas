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

use Definitions, only: iwp

implicit none
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iUHF, iSym, i, nB, nB2, ip, ipO, ipV
character(len=21), parameter :: SecNam = 'RPA_RdOrb_FromRunfile'
integer(kind=iwp), external :: RPA_iUHF

! Restricted (1) or unrestricted (2)
iUHF = RPA_iUHF()

! Allocate memory for CMO
l_CMO(1) = nBas(1)*nOrb(1)
nB2 = nBas(1)**2
do iSym=2,nSym
  l_CMO(1) = l_CMO(1)+nBas(iSym)*nOrb(iSym)
  nB2 = nB2+nBas(iSym)**2
end do
call GetMem('CMO(RPA)','Allo','Real',ip_CMO(1),l_CMO(1))
if (iUHF == 2) then
  l_CMO(2) = l_CMO(1)
  call GetMem('CMO(RPA)','Allo','Real',ip_CMO(2),l_CMO(2))
else
  ip_CMO(2) = 0
  l_CMO(2) = 0
end if

! Read CMO array(s) from Runfile
call Get_CMO(Work(ip_CMO(1)),nB2)
if (iUHF == 2) then
  call Get_dArray('CMO_ab',Work(ip_CMO(2)),nB2)
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
  call GetMem('OccEn','Allo','Real',ip_OccEn(i),l_OccEn(i))
  call GetMem('VirEn','Allo','Real',ip_VirEn(i),l_VirEn(i))
end do
if (iUHF == 1) then
  ip_OccEn(2) = 0
  l_OccEn(2) = 0
  ip_VirEn(2) = 0
  l_VirEn(2) = 0
end if

! Read orbital energies from Runfile
call Get_OrbE(ip_EMO(1),l_EMO(1))
if (l_EMO(1) /= nB) then
  call RPA_Warn(3,SecNam//': unexpected EMO dimension')
end if
ip = ip_EMO(1)
ipO = ip_OccEn(1)
ipV = ip_VirEn(1)
do iSym=1,nSym
  call dCopy_(nOcc(iSym,1),Work(ip),1,Work(ipO),1)
  call dCopy_(nVir(iSym,1),Work(ip+nOcc(iSym,1)),1,Work(ipV),1)
  ip = ip+nOrb(iSym)
  ipO = ipO+nOcc(iSym,1)
  ipV = ipV+nVir(iSym,1)
end do
if (iUHF == 2) then
  l_EMO(2) = l_EMO(1)
  call GetMem('EMO(RPA)','Allo','Real',ip_EMO(2),l_EMO(2))
  call Get_dArray('OrbE_ab',Work(ip_EMO(2)),l_EMO(2))
  ip = ip_EMO(2)
  ipO = ip_OccEn(2)
  ipV = ip_VirEn(2)
  do iSym=1,nSym
    call dCopy_(nOcc(iSym,2),Work(ip),1,Work(ipO),1)
    call dCopy_(nVir(iSym,2),Work(ip+nOcc(iSym,2)),1,Work(ipV),1)
    ip = ip+nOrb(iSym)
    ipO = ipO+nOcc(iSym,2)
    ipV = ipV+nVir(iSym,2)
  end do
end if

end subroutine RPA_RdOrb_FromRunfile
