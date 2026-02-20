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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
!*****************************************************************

subroutine CalcGradCMS(Grad,DDg,lRoots,nSPair)

use Index_Functions, only: iTri
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lRoots, nSPair
real(kind=wp), intent(out) :: Grad(nSPair)
real(kind=wp), intent(in) :: DDg(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: iKL, K, L

do K=2,lRoots
  do L=1,K-1
    iKL = iTri(K-1,L)
    Grad(iKL) = DDg(K,K,K,L)-DDg(L,L,K,L)
  end do
end do
Grad(:) = Two*Grad(:)

return

end subroutine CalcGradCMS
