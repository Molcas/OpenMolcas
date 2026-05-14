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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine Get_D1A(CMO,D1A_MO,D1A_AO,nsym,nbas,nish,nash,ndens)

use Index_Functions, only: iTri
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*), D1A_MO(*)
real(kind=wp), intent(_OUT_) :: D1A_AO(*)
integer(kind=iwp), intent(in) :: nSym, nbas(nsym), nish(nsym), nash(nsym), nDens
integer(kind=iwp) :: i, iAsh, iBas, ii, iIsh, iOff2, iOff3, iSym, j
real(kind=wp), allocatable :: Scr1(:), Tmp1(:,:), Tmp2(:,:)

iOff2 = 1
iOff3 = 1
ii = 0
call mma_allocate(Scr1,2*nDens,Label='Scr1')
do iSym=1,nSym
  iBas = nBas(iSym)
  iAsh = nAsh(iSym)
  iIsh = nIsh(iSym)
  Scr1(iOff3:iOff3+iBas**2-1) = Zero
  if (iAsh /= 0) then
    call mma_allocate(Tmp1,iAsh,iAsh,Label='Tmp1')
    call mma_allocate(Tmp2,iBas,iAsh,Label='Tmp2')
    do i=1,iAsh
      do j=1,iAsh
        Tmp1(j,i) = D1A_MO(iTri(i+ii,j+ii))
      end do
    end do
    ii = ii+iAsh
    call DGEMM_('N','T',iBas,iAsh,iAsh,One,CMO(iOff2+iIsh*iBas),iBas,Tmp1,iAsh,Zero,Tmp2,iBas)
    call DGEMM_('N','T',iBas,iBas,iAsh,One,Tmp2,iBas,CMO(iOff2+iIsh*iBas),iBas,Zero,Scr1(iOff3),iBas)
    call mma_deallocate(Tmp2)
    call mma_deallocate(Tmp1)
  end if
  iOff2 = iOff2+iBas*iBas
  iOff3 = iOff3+iBas*iBas
end do
call Fold2(nSym,nBas,Scr1,D1A_AO)
call mma_deallocate(Scr1)

return

end subroutine Get_D1A
