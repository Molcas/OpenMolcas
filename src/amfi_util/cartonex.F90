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

subroutine cartoneX(L,Lmax,onecontr,ncontrac,MxcontL,onecartX)
!bs arranges the cartesian one-electron-integrals for X  on a quadratic matrix

use Constants, only: Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: L, Lmax, ncontrac, MxcontL
real(kind=wp) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3), onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))
integer(kind=iwp) :: iaddr, icont, jcont, M, Mprime
real(kind=wp) :: pre
!Statement function
integer(kind=iwp) :: ipnt, I, J
ipnt(I,J) = (max(i,j)*(max(i,j)-1))/2+min(i,j)

!bs - + Integrals    m || mprime     mprime=m+1
do Mprime=2,L
  M = mprime-1
  iaddr = ipnt(Mprime+L+1,-M+L+1)
  do jcont=1,ncontrac
    do icont=1,ncontrac
      onecartX(icont,jcont,iaddr) = onecartX(icont,jcont,iaddr)- &
                                    Quart*(onecontr(icont,jcont,Mprime,1)+onecontr(icont,jcont,-Mprime,3))
    end do
  end do
end do
!bs - + Integrals    m || mprime     mprime=m-1
do Mprime=1,L-1
  M = mprime+1
  iaddr = ipnt(Mprime+L+1,-M+L+1)
  do jcont=1,ncontrac
    do icont=1,ncontrac
      onecartX(icont,jcont,iaddr) = onecartX(icont,jcont,iaddr)- &
                                    Quart*(onecontr(icont,jcont,Mprime,3)+onecontr(icont,jcont,-Mprime,1))
    end do
  end do
end do
!bs -1 || 0 integrals
pre = sqrt(0.125_wp)
iaddr = ipnt(L,L+1)
do jcont=1,ncontrac
  do icont=1,ncontrac
    onecartX(icont,jcont,iaddr) = onecartX(icont,jcont,iaddr)-pre*(onecontr(icont,jcont,0,3)+onecontr(icont,jcont,0,1))
  end do
end do

return

end subroutine cartoneX
