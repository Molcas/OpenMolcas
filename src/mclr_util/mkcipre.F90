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
use MCLR_Data, only: ERAS, P1, P1Inv, SS
use input_mclr, only: ERASSCF, lRoots
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, j

call mma_allocate(SS,2*lroots,2*lroots,Label='SS')
SS(:,:) = Zero
do I=1,lroots
  do J=1,lroots
    SS(2*i-1,2*j-1) = P1(iTri(i,j))
  end do
end do
do I=1,lroots
  SS(2*i-1,2*i-1) = SS(2*i-1,2*i-1)+ERAS(I)-ERASSCF(1)
  SS(2*i,2*i-1) = -One
  SS(2*i-1,2*i) = -One
end do
SS(2*lroots-1,2*lroots-1) = SS(2*lroots-1,2*lroots-1)+One
call MatInvert(SS,2*lroots)
do I=1,lroots
  do J=1,lroots
    SS(2*i-1,2*j-1) = SS(2*i-1,2*j-1)+P1INV(itri(i,j))
    SS(2*i,2*j) = SS(2*i,2*j)+P1(iTri(i,j))
  end do
end do
do I=1,lroots
  SS(2*i,2*i-1) = SS(2*i,2*i-1)+One
  SS(2*i-1,2*i) = SS(2*i-1,2*i)+One
end do
call MatInvert(SS,2*lroots)
SS(:,:) = -SS(:,:)
SS(2*lroots,2*lroots) = SS(2*lroots,2*lroots)-One

end subroutine mkcipre
