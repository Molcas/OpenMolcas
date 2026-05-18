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
use constants, only: Zero, One
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp), intent(in) :: nXmo, nSmn, nSaa
real(kind=wp), intent(in) :: Smn(nSmn), Xmo(nXmo)
real(kind=wp), intent(out) :: Saa(nSaa)
real(kind=wp), allocatable :: Z(:)
integer(kind=iwp) mOb, iSym, iX, kX, lX, nBX, j, jK, jX, jZ, lk
real(kind=wp), external :: DDot_

mOb = nBas(1)*nOrb(1)
do iSym=2,nSym
  mOb = max(mOb,nBas(iSym)*nOrb(iSym))
end do
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
