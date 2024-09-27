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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine ConMax(A,nPrim,mPrim,B,nCont,C,mCont)
!***********************************************************************
!                                                                      *
! Object: to find the largest element in the contraction matrix  for   *
!         each primitive index.                                        *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '91                                                 *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrim, mPrim, nCont, mCont
real(kind=wp), intent(out) :: A(nPrim,mPrim)
real(kind=wp), intent(in) :: B(nPrim,nCont), C(mPrim,mCont)
integer(kind=iwp) :: iPrim, jPrim
real(kind=wp) :: Temp
real(kind=wp), external :: DDot_

do iPrim=1,nPrim
  Temp = DDot_(nCont,B(iPrim,1),nPrim,B(iPrim,1),nPrim)
  A(iPrim,:) = Temp
end do

do jPrim=1,mPrim
  Temp = DDot_(mCont,C(jPrim,1),mPrim,C(jPrim,1),mPrim)
  A(:,jPrim) = sqrt(A(:,jPrim)*Temp)
end do

return

end subroutine ConMax
