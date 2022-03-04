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

subroutine cartoneZ(L,Lmax,onecontr,ncontrac,MxcontL,onecartZ)
!bs arranges the cartesian one-electron integrals for Z on a quadratic matrix

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: L, Lmax, ncontrac, MxcontL
real(kind=wp) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3), onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))
integer(kind=iwp) :: iaddr, icont, jcont, Mprime
!Statement function
integer(kind=iwp) :: ipnt, I, J
ipnt(I,J) = (max(i,j)*(max(i,j)-1))/2+min(i,j)

!bs - + Integrals    m || mprime     mprime=m
do Mprime=1,L
  iaddr = ipnt(Mprime+L+1,-mprime+L+1)
  do jcont=1,ncontrac
    do icont=1,ncontrac
      onecartZ(icont,jcont,iaddr) = onecartZ(icont,jcont,iaddr)+ &
                                    Half*(onecontr(icont,jcont,Mprime,2)-onecontr(icont,jcont,-Mprime,2))
    end do
  end do
end do

return

end subroutine cartoneZ
