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

subroutine cartoneY(L,Lmax,onecontr,ncontrac,MxcontL,onecartY)
!bs arranges the cartesian one-electron integrals for Y on a square matrix

use index_functions, only: iTri
use Constants, only: Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, Lmax, ncontrac, MxcontL
real(kind=wp), intent(in) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3)
real(kind=wp), intent(inout) :: onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))
integer(kind=iwp) :: iaddr, M, Mprime
real(kind=wp) :: pre

!bs + + Integrals    m || mprime     mprime=m+1
do Mprime=2,L
  M = mprime-1
  iaddr = iTri(Mprime+L+1,M+L+1)
  onecartY(1:ncontrac,1:ncontrac,iaddr) = onecartY(1:ncontrac,1:ncontrac,iaddr)- &
                                          Quart*(onecontr(1:ncontrac,1:ncontrac,Mprime,1)+onecontr(1:ncontrac,1:ncontrac,-Mprime,3))
end do
!bs - - Integrals    m || mprime     mprime=m-1
do Mprime=1,L-1
  M = mprime+1
  iaddr = iTri(-Mprime+L+1,-M+L+1)
  onecartY(1:ncontrac,1:ncontrac,iaddr) = onecartY(1:ncontrac,1:ncontrac,iaddr)+ &
                                          Quart*(onecontr(1:ncontrac,1:ncontrac,Mprime,3)+onecontr(1:ncontrac,1:ncontrac,-Mprime,1))
end do
!bs  0 || 1 integrals
pre = -sqrt(0.125_wp)
iaddr = iTri(L+1,L+2)
onecartY(1:ncontrac,1:ncontrac,iaddr) = onecartY(1:ncontrac,1:ncontrac,iaddr)+ &
                                        pre*(onecontr(1:ncontrac,1:ncontrac,1,1)+onecontr(1:ncontrac,1:ncontrac,-1,3))

return

end subroutine cartoneY
