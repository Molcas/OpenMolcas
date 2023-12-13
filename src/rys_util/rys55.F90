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

subroutine Rys55(Arg,nArg,Root,Weight,iPntr,nPntr,x0,nMax,R6,R5,R4,R3,R2,R1,R0,W6,W5,W4,W3,W2,W1,W0,ddx,HerW,HerR2,TMax)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             September '90                                            *
!***********************************************************************

use Constants, only: One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, nPntr, iPntr(nPntr), nMax
real(kind=wp), intent(in) :: Arg(nArg), x0(nMax), R6(nMax,5), R5(nMax,5), R4(nMax,5), R3(nMax,5), R2(nMax,5), R1(nMax,5), &
                             R0(nMax,5), W6(nMax,5), W5(nMax,5), W4(nMax,5), W3(nMax,5), W2(nMax,5), W1(nMax,5), W0(nMax,5), ddx, &
                             HerW(5), HerR2(5), TMax
real(kind=wp), intent(out) :: Root(5,nArg), Weight(5,nArg)
integer(kind=iwp) :: iArg, n
real(kind=wp) :: ai, dddx, si, xdInv, z

xdInv = One/ddx
dddx = ddx/Ten+ddx
do iArg=1,nArg
  if (ARg(iArg) < TMax) then
    n = iPntr(int((Arg(iArg)+dddx)*xdInv))
    z = Arg(iArg)-x0(n)
    Root(:,iArg) = (((((R6(n,:)*z+R5(n,:))*z+R4(n,:))*z+R3(n,:))*z+R2(n,:))*z+R1(n,:))*z+R0(n,:)
    Weight(:,iArg) = (((((W6(n,:)*z+W5(n,:))*z+W4(n,:))*z+W3(n,:))*z+W2(n,:))*z+W1(n,:))*z+W0(n,:)
  else
    ai = One/Arg(iArg)
    si = sqrt(ai)
    Root(:,iArg) = HerR2(:)*ai
    Weight(:,iArg) = HerW(:)*si
  end if
end do

return

end subroutine Rys55
