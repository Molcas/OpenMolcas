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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Header
logical(kind=iwp), intent(in) :: PrOcc, PrEne
real(kind=wp), intent(in) :: ThrOcc, ThrEne, Ene(*), Occ(*), CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), iPrForm, IndxT(*)
character(len=*), intent(in) :: Nme(*)
integer(kind=iwp) :: ik, iSym, k, k1, k2, kk, l_CMO, nOrb_(8), nOrbT
real(kind=wp), allocatable :: CMO_(:), EOr_(:), Occ_(:)

nOrb_(1:nSym) = nOrb(:)
kk = 0
do iSym=1,nSym
  do k=1,nBas(iSym)
    ik = kk+k
    if (IndxT(ik) == 7) nOrb_(iSym) = nOrb_(iSym)-1
  end do
  kk = kk+nBas(iSym)
end do

nOrbT = 0
do iSym=1,nSym
  nOrbT = nOrbT+nOrb_(iSym)
end do

l_CMO = 0
do iSym=1,nSym
  l_CMO = l_CMO+nBas(iSym)*nOrb_(iSym)
end do

call mma_allocate(CMO_,l_CMO,label='CMO_')
call mma_allocate(Occ_,nOrbT,label='Occ_')
call mma_allocate(EOr_,nOrbT,label='Eor_')

k1 = 1
k2 = 1
do iSym=1,nSym
  call dCopy_(nBas(iSym)*nOrb_(iSym),CMO(k1),1,CMO_(k2),1)
  k1 = k1+nBas(iSym)*nBas(iSym)
  k2 = k2+nBas(iSym)*nOrb_(iSym)
end do

k1 = 1
k2 = 1
do iSym=1,nSym
  call dCopy_(nOrb_(iSym),Occ(k1),1,Occ_(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb_(iSym)
end do

k1 = 1
k2 = 1
do iSym=1,nSym
  call dCopy_(nOrb_(iSym),Ene(k1),1,EOr_(k2),1)
  k1 = k1+nBas(iSym)
  k2 = k2+nOrb_(iSym)
end do

call PriMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb_,Nme,EOr_,Occ_,CMO_,iPrForm)

call mma_deallocate(CMO_)
call mma_deallocate(Occ_)
call mma_deallocate(EOr_)

end subroutine PriMO_Localisation
