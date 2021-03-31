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

subroutine WrVec_Localisation(FName,Lu,Label,nSym,nBas,nOrb,CMO,Occ,EOrb,IndT,Title)
! Thomas Bondo Pedersen, July 2010.
!
! Write orbital info.
! This is a work-around to fix bugs when orbitals are deleted.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=6), intent(in) :: FName
integer(kind=iwp), intent(in) :: Lu, nSym, nBas(nSym), nOrb(nSym), IndT(*)
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: CMO(*), Occ(*), EOrb(*)
character(len=*), intent(inout) :: Title
#include "WrkSpc.fh"
integer(kind=iwp) :: ip_CMO, ip_EOr, ip_Ind, ip_Occ, iSym, k1, k2, l_CMO, l_EOr, l_Ind, l_Occ
logical(kind=iwp) :: Write_CMO, Write_EOr, Write_Ind, Write_Occ

Write_CMO = index(Label,'C') /= 0
Write_Occ = index(Label,'O') /= 0
Write_EOr = index(Label,'E') /= 0
Write_Ind = index(Label,'I') /= 0

if (Write_CMO) then
  l_CMO = nBas(1)*nOrb(1)
  do iSym=2,nSym
    l_CMO = l_CMO+nBas(iSym)*nOrb(iSym)
  end do
  call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
  k1 = 1
  k2 = ip_CMO
  do iSym=1,nSym
    call dCopy_(nBas(iSym)*nOrb(iSym),CMO(k1),1,Work(k2),1)
    k1 = k1+nBas(iSym)*nBas(iSym)
    k2 = k2+nBas(iSym)*nOrb(iSym)
  end do
else
  l_CMO = 1
  call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
  Work(ip_CMO) = Zero
end if

if (Write_Occ) then
  l_Occ = nOrb(1)
  do iSym=2,nSym
    l_Occ = l_Occ+nOrb(iSym)
  end do
  call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
  k1 = 1
  k2 = ip_Occ
  do iSym=1,nSym
    call dCopy_(nOrb(iSym),Occ(k1),1,Work(k2),1)
    k1 = k1+nBas(iSym)
    k2 = k2+nOrb(iSym)
  end do
else
  l_Occ = 1
  call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
  Work(ip_Occ) = Zero
end if

if (Write_EOr) then
  l_EOr = nOrb(1)
  do iSym=2,nSym
    l_EOr = l_EOr+nOrb(iSym)
  end do
  call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
  k1 = 1
  k2 = ip_EOr
  do iSym=1,nSym
    call dCopy_(nOrb(iSym),EOrb(k1),1,Work(k2),1)
    k1 = k1+nBas(iSym)
    k2 = k2+nOrb(iSym)
  end do
else
  l_EOr = 1
  call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
  Work(ip_EOr) = Zero
end if

if (Write_Ind) then
  l_Ind = 56
  call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)
  call iCopy(l_Ind,IndT,1,iWork(ip_Ind),1)
else
  l_Ind = 1
  call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)
  iWork(ip_Ind) = 0
end if

call WrVec(FName,Lu,Label,nSym,nBas,nOrb,Work(ip_CMO),Work(ip_Occ),Work(ip_EOr),iWork(ip_Ind),Title)

call GetMem('Ind_','Free','Inte',ip_Ind,l_Ind)
call GetMem('EOr_','Free','Real',ip_EOr,l_EOr)
call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

end subroutine WrVec_Localisation
