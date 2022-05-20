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

subroutine Rys66(Arg,nArg,Root,Weight,iPntr,nPntr,x0,nMax,R6,R5,R4,R3,R2,R1,R0,W6,W5,W4,W3,W2,W1,W0,ddx,HerW,HerR2,TMax)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             September '90                                            *
!***********************************************************************

use Constants, only: One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nArg, nPntr, iPntr(nPntr), nMax
real(kind=wp) :: Arg(nArg), Root(6,nArg), Weight(6,nArg), x0(nMax), R6(nMax,6), R5(nMax,6), R4(nMax,6), R3(nMax,6), R2(nMax,6), &
                 R1(nMax,6), R0(nMax,6), W6(nMax,6), W5(nMax,6), W4(nMax,6), W3(nMax,6), W2(nMax,6), W1(nMax,6), W0(nMax,6), ddx, &
                 HerW(6), HerR2(6), TMax
integer(kind=iwp) :: iArg, n
real(kind=wp) :: ai, dddx, si, xdInv, z

xdInv = One/ddx
dddx = ddx/Ten+ddx
do iArg=1,nArg
  if (Arg(iArg) < TMax) then
    n = iPntr(int((Arg(iArg)+dddx)*xdInv))
    z = Arg(iArg)-x0(n)
    Root(1,iArg) = (((((R6(n,1)*z+R5(n,1))*z+R4(n,1))*z+R3(n,1))*z+R2(n,1))*z+R1(n,1))*z+R0(n,1)
    Root(2,iArg) = (((((R6(n,2)*z+R5(n,2))*z+R4(n,2))*z+R3(n,2))*z+R2(n,2))*z+R1(n,2))*z+R0(n,2)
    Root(3,iArg) = (((((R6(n,3)*z+R5(n,3))*z+R4(n,3))*z+R3(n,3))*z+R2(n,3))*z+R1(n,3))*z+R0(n,3)
    Root(4,iArg) = (((((R6(n,4)*z+R5(n,4))*z+R4(n,4))*z+R3(n,4))*z+R2(n,4))*z+R1(n,4))*z+R0(n,4)
    Root(5,iArg) = (((((R6(n,5)*z+R5(n,5))*z+R4(n,5))*z+R3(n,5))*z+R2(n,5))*z+R1(n,5))*z+R0(n,5)
    Root(6,iArg) = (((((R6(n,6)*z+R5(n,6))*z+R4(n,6))*z+R3(n,6))*z+R2(n,6))*z+R1(n,6))*z+R0(n,6)
    Weight(1,iArg) = (((((W6(n,1)*z+W5(n,1))*z+W4(n,1))*z+W3(n,1))*z+W2(n,1))*z+W1(n,1))*z+W0(n,1)
    Weight(2,iArg) = (((((W6(n,2)*z+W5(n,2))*z+W4(n,2))*z+W3(n,2))*z+W2(n,2))*z+W1(n,2))*z+W0(n,2)
    Weight(3,iArg) = (((((W6(n,3)*z+W5(n,3))*z+W4(n,3))*z+W3(n,3))*z+W2(n,3))*z+W1(n,3))*z+W0(n,3)
    Weight(4,iArg) = (((((W6(n,4)*z+W5(n,4))*z+W4(n,4))*z+W3(n,4))*z+W2(n,4))*z+W1(n,4))*z+W0(n,4)
    Weight(5,iArg) = (((((W6(n,5)*z+W5(n,5))*z+W4(n,5))*z+W3(n,5))*z+W2(n,5))*z+W1(n,5))*z+W0(n,5)
    Weight(6,iArg) = (((((W6(n,6)*z+W5(n,6))*z+W4(n,6))*z+W3(n,6))*z+W2(n,6))*z+W1(n,6))*z+W0(n,6)
  else
    ai = One/Arg(iArg)
    si = sqrt(ai)
    Root(1,iArg) = HerR2(1)*ai
    Root(2,iArg) = HerR2(2)*ai
    Root(3,iArg) = HerR2(3)*ai
    Root(4,iArg) = HerR2(4)*ai
    Root(5,iArg) = HerR2(5)*ai
    Root(6,iArg) = HerR2(6)*ai
    Weight(1,iArg) = HerW(1)*si
    Weight(2,iArg) = HerW(2)*si
    Weight(3,iArg) = HerW(3)*si
    Weight(4,iArg) = HerW(4)*si
    Weight(5,iArg) = HerW(5)*si
    Weight(6,iArg) = HerW(6)*si
  end if
end do

return

end subroutine Rys66
