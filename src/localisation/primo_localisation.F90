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

implicit none
#include "Molcas.fh"
character*(*) Header
character*(LENIN8) Nme(*)
logical PrOcc, PrEne
real*8 ThrOcc, ThrEne
integer nSym
integer nBas(nSym), nOrb(nSym)
real*8 Ene(*), Occ(*), CMO(*)
integer iPrForm
integer IndxT(*)
#include "WrkSpc.fh"

integer ip_CMO, l_CMO
integer ip_Occ, l_OCc
integer ip_EOr, l_EOr
integer nOrbT, iSym, k1, k2
integer k, kk, ik, nOrb_(8)

call Icopy(nSym,nOrb,1,nOrb_,1)
kk = 0
do iSym=1,nSym
  do k=1,nBas(iSym)
    ik = kk+k
    if (IndxT(ik) == 7) nOrb(iSym) = nOrb(iSym)-1
  end do
  kk = kk+nBas(iSym)
end do

nOrbT = nOrb(1)
do iSym=2,nSym
  nOrbT = nOrbT+nOrb(iSym)
end do

l_CMO = nBas(1)*nOrb(1)
do iSym=2,nSym
  l_CMO = l_CMO+nBas(iSym)*nOrb(iSym)
end do
l_Occ = nOrbT
l_EOr = nOrbT

call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
call GetMem('Eor_','Allo','Real',ip_EOr,l_EOr)

k1 = 1
k2 = ip_CMO
do iSym=1,nSym
  call dCopy_(nBas(iSym)*nOrb(iSym),CMO(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)*nBas(iSym)
  k2 = k2+nBas(iSym)*nOrb(iSym)
end do

k1 = 1
k2 = ip_Occ
do iSym=1,nSym
  call dCopy_(nOrb(iSym),Occ(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb(iSym)
end do

k1 = 1
k2 = ip_EOr
do iSym=1,nSym
  call dCopy_(nOrb(iSym),Ene(k1),1,Work(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb(iSym)
end do

call PriMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,Nme,Work(ip_EOr),Work(ip_Occ),Work(ip_CMO),iPrForm)

call GetMem('Eor_','Free','Real',ip_EOr,l_EOr)
call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

call Icopy(nSym,nOrb_,1,nOrb,1)

end subroutine PriMO_Localisation
