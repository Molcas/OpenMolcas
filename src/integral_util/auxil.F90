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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Auxil(T,nT,Fm,mHigh)
!***********************************************************************
!                                                                      *
!     Object: to compute the auxiliary functions in quadruple precision*
!             for a number of arguments.                               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, mHigh
real(kind=wp), intent(in) :: T(nT)
real(kind=wp), intent(out) :: Fm(nT,0:mHigh)
integer(kind=iwp) :: i, m
real(kind=wp) :: Ti

call HighFm(Fm(:,mHigh),T,mHigh,nT)

! Now use recursion formula for Fm, 0<=m<mHigh

do i=1,nT
  Ti = T(i)
  do m=mHigh-1,0,-1
    Fm(i,m) = (Two*Ti*Fm(i,m+1)+exp(-Ti))/real(2*m+1,kind=wp)
  end do
end do
!call RecPrt(' Fm',' ',Fm,nT,mHigh+1)

return

end subroutine Auxil
