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
!               2011, Francesco Aquilante                              *
!***********************************************************************

subroutine PriMO_Localisation(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,Nme,Ene,Occ,CMO,iPrForm,IndxT)
! Thomas Bondo Pedersen, July 2010.
!
! Print MOs.
! This is a work-around to fix bugs when orbitals are deleted.
!
! F. Aquilante, Nov 2011   (Print only non-deleted orbitals)

use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
character(len=*), intent(in) :: Header
logical(kind=iwp), intent(in) :: PrOcc, PrEne
real(kind=wp), intent(in) :: ThrOcc, ThrEne, Ene(*), Occ(*), CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), iPrForm, IndxT(*)
character(len=LenIn8), intent(in) :: Nme(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: ik, ip_CMO, ip_EOr, ip_Occ, iSym, k, k1, k2, kk, l_CMO, l_EOr, l_Occ, nOrb_(8), nOrbT

call Icopy(nSym,nOrb,1,nOrb_,1)
kk = 0
do iSym=1,nSym
  do k=1,nBas(iSym)
    ik = kk+k
    if (IndxT(ik) == 7) nOrb_(iSym) = nOrb_(iSym)-1
  end do
  kk = kk+nBas(iSym)
end do

nOrbT = nOrb_(1)
do iSym=2,nSym
  nOrbT = nOrbT+nOrb_(iSym)
end do

l_CMO = nBas(1)*nOrb_(1)
do iSym=2,nSym
  l_CMO = l_CMO+nBas(iSym)*nOrb_(iSym)
end do
l_Occ = nOrbT
l_EOr = nOrbT

call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
call GetMem('Eor_','Allo','Real',ip_EOr,l_EOr)

k1 = 1
k2 = ip_CMO
do iSym=1,nSym
  call dCopy_(nBas(iSym)*nOrb_(iSym),CMO(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)*nBas(iSym)
  k2 = k2+nBas(iSym)*nOrb_(iSym)
end do

k1 = 1
k2 = ip_Occ
do iSym=1,nSym
  call dCopy_(nOrb_(iSym),Occ(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb_(iSym)
end do

k1 = 1
k2 = ip_EOr
do iSym=1,nSym
  call dCopy_(nOrb_(iSym),Ene(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb_(iSym)
end do

call PriMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb_,Nme,Work(ip_EOr),Work(ip_Occ),Work(ip_CMO),iPrForm)

call GetMem('Eor_','Free','Real',ip_EOr,l_EOr)
call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

end subroutine PriMO_Localisation
