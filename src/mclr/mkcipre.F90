!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine mkcipre()

use Index_Functions, only: iTri
use MCLR_Data, only: SS, ERAS, P1, P1Inv
use input_mclr, only: lRoots, ERASSCF
use stdalloc, only: mma_allocate
use Constants, only: One

implicit none
integer i, j, iRec
! Statement functions
irec(i,j) = i+(j-1)*2*lroots

call mma_allocate(SS,4*lroots**2,Label='SS')
do I=1,lroots
  do J=1,lroots
    SS(irec(2*i-1,2*j-1)) = P1(iTri(i,j))
  end do
end do
do I=1,lroots
  SS(irec(2*i-1,2*i-1)) = SS(irec(2*i-1,2*i-1))+ERAS(I)-ERASSCF(1)
  SS(irec(2*i,2*i-1)) = -One
  SS(irec(2*i-1,2*i)) = -One
end do
SS(irec(2*lroots-1,2*lroots-1)) = SS(irec(2*lroots-1,2*lroots-1))+One
call MatInvert(SS,2*lroots)
do I=1,lroots
  do J=1,lroots
    SS(irec(2*i-1,2*j-1)) = SS(irec(2*i-1,2*j-1))+P1INV(itri(i,j))
    SS(irec(2*i,2*j)) = SS(irec(2*i,2*j))+P1(iTri(i,j))
  end do
end do
do I=1,lroots
  SS(irec(2*i,2*i-1)) = SS(irec(2*i,2*i-1))+One
  SS(irec(2*i-1,2*i)) = SS(irec(2*i-1,2*i))+One
end do
call MatInvert(SS,2*lroots)
call DSCAL_(4*lroots**2,-One,SS,1)
SS(irec(2*lroots,2*lroots)) = SS(irec(2*lroots,2*lroots))-One

end subroutine mkcipre
