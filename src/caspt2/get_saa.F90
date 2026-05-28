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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine get_Saa(nSym,nBas,nOrb,Smn,nSmn,Xmo,nXmo,Saa,nSaa)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), nSmn, nXmo, nSaa
real(kind=wp), intent(in) :: Smn(nSmn), Xmo(nXmo)
real(kind=wp), intent(out) :: Saa(nSaa)
integer(kind=iwp) :: iSym, iX, j, jK, jX, jZ, kX, lk, lX, mOb, nBX
real(kind=wp), allocatable :: Z(:)
real(kind=wp), external :: DDot_

mOb = maxval(nBas(:)*nOrb(:))
call mma_allocate(Z,mOb,Label='Z')

iX = 1
kX = 1
lX = 1
do iSym=1,nSym
  nBx = max(1,nBas(iSym))
  call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),One,Smn(iX),nBx,Xmo(kX),nBx,Zero,Z,nBx)
  do j=0,nOrb(iSym)-1
    jK = nBas(iSym)*j
    lk = kX+jK
    jZ = 1+jK
    jX = lX+j
    Saa(jX) = ddot_(nBas(iSym),Xmo(lk),1,Z(jZ),1)
  end do
  iX = iX+nBas(iSym)**2
  kX = kX+nBas(iSym)*nOrb(iSym)
  lX = lX+nOrb(iSym)
end do
call mma_deallocate(Z)

end subroutine get_Saa
