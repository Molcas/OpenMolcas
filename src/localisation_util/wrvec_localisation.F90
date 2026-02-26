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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=6), intent(in) :: FName
integer(kind=iwp), intent(in) :: Lu, nSym, nBas(nSym), nOrb(nSym), IndT(*)
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: CMO(*), Occ(*), EOrb(*)
character(len=*), intent(inout) :: Title
integer(kind=iwp) :: iSym, k1, k2, l_CMO, l_EOr, l_Ind, l_Occ
logical(kind=iwp) :: Write_CMO, Write_EOr, Write_Ind, Write_Occ
integer(kind=iwp), allocatable :: Ind_(:)
real(kind=wp), allocatable :: CMO_(:), Occ_(:), EOr_(:)

Write_CMO = index(Label,'C') /= 0
Write_Occ = index(Label,'O') /= 0
Write_EOr = index(Label,'E') /= 0
Write_Ind = index(Label,'I') /= 0

if (Write_CMO) then
  l_CMO = nBas(1)*nOrb(1)
  do iSym=2,nSym
    l_CMO = l_CMO+nBas(iSym)*nOrb(iSym)
  end do
  call mma_allocate(CMO_,l_CMO,label='CMO')
  k1 = 1
  k2 = 1
  do iSym=1,nSym
    call dCopy_(nBas(iSym)*nOrb(iSym),CMO(k1),1,CMO_(k2),1)
    k1 = k1+nBas(iSym)*nBas(iSym)
    k2 = k2+nBas(iSym)*nOrb(iSym)
  end do
else
  call mma_allocate(CMO_,1,label='CMO')
  CMO_(1) = Zero
end if

if (Write_Occ) then
  l_Occ = nOrb(1)
  do iSym=2,nSym
    l_Occ = l_Occ+nOrb(iSym)
  end do
  call mma_allocate(Occ_,l_Occ,label='Occ_')
  k1 = 1
  k2 = 1
  do iSym=1,nSym
    call dCopy_(nOrb(iSym),Occ(k1),1,Occ_(k2),1)
    k1 = k1+nBas(iSym)
    k2 = k2+nOrb(iSym)
  end do
else
  call mma_allocate(Occ_,1,label='Occ_')
  Occ_(1) = Zero
end if

if (Write_EOr) then
  l_EOr = nOrb(1)
  do iSym=2,nSym
    l_EOr = l_EOr+nOrb(iSym)
  end do
  call mma_allocate(EOr_,l_EOr,label='EOr')
  k1 = 1
  k2 = 1
  do iSym=1,nSym
    call dCopy_(nOrb(iSym),EOrb(k1),1,EOr_(k2),1)
    k1 = k1+nBas(iSym)
    k2 = k2+nOrb(iSym)
  end do
else
  call mma_allocate(EOr_,1,label='EOr')
  EOr_(1) = Zero
end if

if (Write_Ind) then
  l_Ind = 7*8
  call mma_allocate(Ind_,l_Ind,label='Ind_')
  Ind_(:) = IndT(1:l_Ind)
else
  call mma_allocate(Ind_,1,label='Ind_')
  Ind_(1) = 0
end if

call WrVec(FName,Lu,Label,nSym,nBas,nOrb,CMO_,Occ_,EOr_,Ind_,Title)

call mma_deallocate(CMO_)
call mma_deallocate(Occ_)
call mma_deallocate(EOr_)
call mma_deallocate(Ind_)

end subroutine WrVec_Localisation
